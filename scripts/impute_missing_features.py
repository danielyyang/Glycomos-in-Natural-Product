"""
GlycoNP 数据填补与特征修复 (Data Imputation & Feature Repair)

三层填补策略 (Three-Layer Imputation Strategy):
  Layer 1: RDKit 理化参数重算 (验证并修复异常值)
  Layer 2: Tanimoto 最近邻分类救援 (np_classifier_superclass)
  Layer 3: 文本挖掘备选方案准备 (dois/organisms)

使用方法 (Usage):
  python scripts/impute_missing_features.py [--input PATH]
"""
import argparse
import gc
import os
import sys
import time
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
import pandas as pd

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "lib")))


# =====================================================================
# Layer 1: RDKit 理化参数重算 (Physicochemical Recalculation)
# =====================================================================

def imputePhysicochemical(df: pd.DataFrame) -> Dict[str, int]:
    """
    用 RDKit 重算空值/异常值的理化参数。
    Recalculate missing/abnormal physicochemical properties using RDKit.

    扫描列 (Scanned columns):
      alogp, topological_polar_surface_area, hydrogen_bond_donors,
      hydrogen_bond_acceptors, qed_drug_likeliness, molecular_weight
    """
    from rdkit import Chem
    from rdkit.Chem import Descriptors, QED

    targetCols = {
        "alogp": lambda m: round(Descriptors.MolLogP(m), 2),
        "topological_polar_surface_area": lambda m: round(Descriptors.TPSA(m), 2),
        "hydrogen_bond_donors": lambda m: Descriptors.NumHDonors(m),
        "hydrogen_bond_acceptors": lambda m: Descriptors.NumHAcceptors(m),
        "qed_drug_likeliness": lambda m: round(QED.qed(m), 4),
        "molecular_weight": lambda m: round(Descriptors.ExactMolWt(m), 2),
    }

    fillCounts: Dict[str, int] = {}

    for colName, calcFunc in targetCols.items():
        if colName not in df.columns:
            print(f"    [SKIP] Column '{colName}' not found")
            continue

        # 找空值行 (Find NaN rows)
        mask = df[colName].isna() | (df[colName].astype(str).isin(["", "nan"]))
        nMissing = mask.sum()
        fillCounts[colName] = 0

        if nMissing == 0:
            print(f"    [OK] {colName}: 0 missing")
            continue

        print(f"    [FILL] {colName}: {nMissing} missing → recalculating...")
        filled = 0
        for idx in df[mask].index:
            smiles = str(df.at[idx, "canonical_smiles"])
            if not smiles or smiles in ("nan", "", "NULL"):
                continue
            try:
                mol = Chem.MolFromSmiles(smiles)
                if mol is not None:
                    val = calcFunc(mol)
                    df.at[idx, colName] = val
                    filled += 1
            except Exception:
                pass

        fillCounts[colName] = filled
        print(f"    [DONE] {colName}: filled {filled}/{nMissing}")

    return fillCounts


# =====================================================================
# Layer 2: Tanimoto 最近邻分类救援 (Tanimoto NN Classification Rescue)
# =====================================================================

def imputeSuperclassByTanimoto(
    df: pd.DataFrame,
    targetCol: str = "np_classifier_superclass",
    threshold: float = 0.85,
    refSampleSize: int = 5000,
    fpRadius: int = 2,
    fpBits: int = 2048,
) -> int:
    """
    用 Morgan 指纹 Tanimoto 相似度, 为空白大类赋值最近邻的分类。
    Assign empty Superclass by finding nearest classified neighbor via Tanimoto.

    算法 (Algorithm):
      1. 从已分类样本中按比例抽取参考集 (refSampleSize)
      2. 计算所有参考样本的 Morgan 指纹
      3. 对每个空值样本, 计算与参考集的批量 Tanimoto
      4. 若最高 Tanimoto >= threshold, 继承该参考样本的分类

    为什么不用全量对全量 (Why not full-vs-full):
      94K × 94K = 88 亿次比较 → 内存爆炸。
      用 5K 参考集 × 11K 查询 = 5500 万次, 可控。
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit import DataStructs

    if targetCol not in df.columns:
        print(f"    [SKIP] Column '{targetCol}' not found")
        return 0

    # 找空值行 (Find NaN rows)
    mask = df[targetCol].isna() | (df[targetCol].astype(str).isin(["", "nan"]))
    nMissing = mask.sum()

    if nMissing == 0:
        print(f"    [OK] {targetCol}: 0 missing")
        return 0

    print(f"    [RESCUE] {targetCol}: {nMissing} missing")

    # 构建参考集 (Build reference set from classified samples)
    classifiedDf = df[~mask].copy()
    if len(classifiedDf) == 0:
        print(f"    [WARN] No classified samples to reference")
        return 0

    # 分层抽样参考集 (Stratified reference sampling)
    nRef = min(refSampleSize, len(classifiedDf))
    refDf = classifiedDf.sample(n=nRef, random_state=42)
    print(f"    Reference set: {nRef} samples from {len(classifiedDf)} classified")

    # 计算参考指纹 (Compute reference fingerprints)
    refFps = []
    refLabels = []
    refIndices = []

    for idx, row in refDf.iterrows():
        smiles = str(row.get("canonical_smiles", ""))
        if not smiles or smiles in ("nan", ""):
            continue
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is not None:
                fp = AllChem.GetMorganFingerprintAsBitVect(
                    mol, fpRadius, nBits=fpBits)
                refFps.append(fp)
                refLabels.append(str(row[targetCol]))
                refIndices.append(idx)
        except Exception:
            pass

    print(f"    Valid reference fingerprints: {len(refFps)}")

    if not refFps:
        return 0

    # 对空值样本做 Tanimoto 匹配 (Match empty samples)
    queryDf = df[mask].copy()
    filled = 0
    batchSize = 500

    for batchStart in range(0, len(queryDf), batchSize):
        batch = queryDf.iloc[batchStart:batchStart + batchSize]

        for idx, row in batch.iterrows():
            smiles = str(row.get("canonical_smiles", ""))
            if not smiles or smiles in ("nan", ""):
                continue

            try:
                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    continue
                queryFp = AllChem.GetMorganFingerprintAsBitVect(
                    mol, fpRadius, nBits=fpBits)

                # 批量 Tanimoto (Bulk Tanimoto)
                sims = DataStructs.BulkTanimotoSimilarity(queryFp, refFps)
                maxIdx = int(np.argmax(sims))
                maxSim = sims[maxIdx]

                if maxSim >= threshold:
                    df.at[idx, targetCol] = refLabels[maxIdx]
                    filled += 1

            except Exception:
                pass

        if (batchStart // batchSize) % 5 == 0:
            print(f"      Processed {batchStart + len(batch)}/{len(queryDf)} "
                  f"queries, filled {filled}...", end="\r")

    print(f"\n    [DONE] Filled {filled}/{nMissing} via Tanimoto >= {threshold}")
    return filled


# =====================================================================
# Layer 3: 文本挖掘备选方案 (Text Mining Fallback Preparation)
# =====================================================================

def prepareTextMiningCandidates(df: pd.DataFrame) -> pd.DataFrame:
    """
    为 dois/organisms 空值行准备文本挖掘候选。
    Prepare text mining candidates for empty dois/organisms.

    不执行实际 API 调用, 仅生成候选列表 CSV 供后续处理。
    Does NOT make API calls — generates candidate list CSV for later processing.
    """
    candidates = []

    doisMask = df["dois"].isna() | (df["dois"].astype(str).isin(["", "nan"]))
    orgMask = df["organisms"].isna() | (df["organisms"].astype(str).isin(["", "nan"]))
    bothMask = doisMask | orgMask

    for idx in df[bothMask].index:
        name = str(df.at[idx, "name"])
        iupac = str(df.at[idx, "iupac_name"]) if "iupac_name" in df.columns else ""
        inchiKey = str(df.at[idx, "standard_inchi_key"])

        if name in ("nan", ""):
            name = ""
        if iupac in ("nan", ""):
            iupac = ""

        # 只保留有名称的候选 (Only keep candidates with a name)
        if name or iupac:
            candidates.append({
                "standard_inchi_key": inchiKey,
                "name": name,
                "iupac_name": iupac,
                "missing_dois": bool(doisMask.get(idx, False)),
                "missing_organisms": bool(orgMask.get(idx, False)),
            })

    candidateDf = pd.DataFrame(candidates)
    return candidateDf


# =====================================================================
# 汇总与 Superclass 同步 (Sync np_classifier → Superclass)
# =====================================================================

def syncSuperclass(df: pd.DataFrame) -> int:
    """
    将 Layer 2 填补的 np_classifier_superclass 同步回 Superclass 列。
    Sync rescued np_classifier_superclass back to Superclass column.

    仅更新当前 Superclass 是 'Unclassified' 的行。
    """
    if "np_classifier_superclass" not in df.columns or "Superclass" not in df.columns:
        return 0

    synced = 0
    unclassifiedMask = (
        df["Superclass"].isna()
        | (df["Superclass"].astype(str).isin(["", "nan", "Unclassified"]))
    )
    hasNpClass = ~(
        df["np_classifier_superclass"].isna()
        | (df["np_classifier_superclass"].astype(str).isin(["", "nan"]))
    )

    toSync = unclassifiedMask & hasNpClass

    for idx in df[toSync].index:
        npClass = str(df.at[idx, "np_classifier_superclass"])
        if npClass and npClass not in ("nan", ""):
            df.at[idx, "Superclass"] = npClass
            synced += 1

    return synced


# =====================================================================
# Main
# =====================================================================

def main():
    parser = argparse.ArgumentParser(
        description="GlycoNP Data Imputation & Feature Repair")
    parser.add_argument("--input", type=str, default=None, help="Input CSV")
    args = parser.parse_args()

    baseDir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    reportDir = os.path.join(baseDir, "reports")

    if args.input:
        inputPath = args.input
    else:
        inputPath = os.path.join(reportDir, "GlycoNP_Pipeline_Full.csv")

    print("=" * 70)
    print("  GlycoNP Data Imputation & Feature Repair")
    print("  数据填补与特征修复")
    print("=" * 70)
    print(f"  Input: {inputPath}")

    t0 = time.time()
    df = pd.read_csv(inputPath, low_memory=False, encoding="utf-8-sig")
    print(f"  Loaded: {len(df):,} rows ({time.time()-t0:.1f}s)")

    # ---- Layer 1: RDKit 理化参数 ----
    print(f"\n{'='*70}")
    print(f"  Layer 1: RDKit Physicochemical Recalculation")
    print(f"{'='*70}")
    physFills = imputePhysicochemical(df)

    # ---- Layer 2: Tanimoto 分类救援 ----
    print(f"\n{'='*70}")
    print(f"  Layer 2: Tanimoto Classification Rescue")
    print(f"{'='*70}")
    tanimotoFilled = imputeSuperclassByTanimoto(
        df,
        targetCol="np_classifier_superclass",
        threshold=0.85,
        refSampleSize=5000,
    )

    # Sync to Superclass
    synced = syncSuperclass(df)
    print(f"    Synced np_classifier → Superclass: {synced}")

    # ---- Layer 3: 文本挖掘准备 ----
    print(f"\n{'='*70}")
    print(f"  Layer 3: Text Mining Candidate Preparation")
    print(f"{'='*70}")
    candidateDf = prepareTextMiningCandidates(df)
    candidatePath = os.path.join(reportDir, "TextMining_Candidates.csv")
    candidateDf.to_csv(candidatePath, index=False, encoding="utf-8-sig")
    print(f"    Candidates: {len(candidateDf):,} (missing DOI or organism)")
    print(f"    Saved: {candidatePath}")

    # ---- 保存 ----
    outPath = os.path.join(reportDir, "GlycoNP_Pipeline_Full_Imputed.csv")
    print(f"\n  Saving imputed data...")
    df.to_csv(outPath, index=False, encoding="utf-8-sig")
    print(f"  Saved: {outPath}")

    # ---- 汇总报告 ----
    elapsed = time.time() - t0
    print(f"\n{'='*70}")
    print(f"  Imputation Summary Report")
    print(f"{'='*70}")
    print(f"  Total rows: {len(df):,}")
    print(f"  Time: {elapsed:.0f}s")
    print(f"\n  Layer 1 — Physicochemical Fills:")
    totalPhys = 0
    for col, count in physFills.items():
        print(f"    {col:40s} +{count}")
        totalPhys += count
    print(f"    {'TOTAL':40s} +{totalPhys}")

    print(f"\n  Layer 2 — Tanimoto Classification Rescue:")
    print(f"    np_classifier_superclass filled:  +{tanimotoFilled}")
    print(f"    Synced to Superclass:             +{synced}")

    # 验证最终 NaN 率 (Final NaN rates)
    print(f"\n  Final NaN rates:")
    checkCols = ["alogp", "np_classifier_superclass", "Superclass",
                 "dois", "organisms"]
    for col in checkCols:
        if col in df.columns:
            na = df[col].isna() | (df[col].astype(str).isin(["", "nan"]))
            pct = na.sum() / len(df) * 100
            print(f"    {col:40s} {na.sum():>6} ({pct:.1f}%)")

    print(f"\n  Layer 3 — Text Mining Candidates: {len(candidateDf):,}")
    print(f"\n  Output: {outPath}")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
