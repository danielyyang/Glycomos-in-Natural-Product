"""
GlycoNP 生物活性与文献概念挖掘 (Bioactivity & Literature Mining)

两大数据引擎 (Two Data Engines):
  Engine 1: ChEMBL 本地静态合并 — InChIKey → Targets + IC50/Ki
  Engine 2: OpenAlex 合规文献概念提取 — DOI → NLP Concepts (Polite Pool)

最终产物: 糖链-活性映射表 (Sugar–Activity Mapping)

使用方法 (Usage):
  python scripts/mine_bioactivity_and_literature.py [--pilot 500] [--email user@example.com]
"""
import argparse
import json
import os
import sys
import time
from collections import Counter
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))


# =====================================================================
# Engine 1: ChEMBL 本地静态合并 (ChEMBL Local Static Merge)
# =====================================================================

def mergeChembl(
    df: pd.DataFrame,
    chemblMolPath: str,
    chemblActPath: str,
) -> pd.DataFrame:
    """
    通过 InChIKey 将 GlycoNP 与 ChEMBL 本地数据进行零网络请求合并。
    Merge GlycoNP with ChEMBL via InChIKey — zero network requests.

    合并策略 (Merge Strategy):
      1. 加载 chembl_molecules.csv → 提取 (standard_inchi_key, molecule_chembl_id)
      2. Inner join GlycoNP.standard_inchi_key ⋈ ChEMBL.standard_inchi_key
      3. 加载 chembl_activities.csv → 仅保留命中 chembl_id 的行
      4. 提取 target_pref_name, standard_type, pchembl_value
    """
    print(f"  [ChEMBL] Loading molecules: {chemblMolPath}")
    molDf = pd.read_csv(chemblMolPath, dtype=str, low_memory=False,
                        usecols=["molecule_chembl_id", "standard_inchi_key"])
    molDf = molDf.dropna(subset=["standard_inchi_key"])
    print(f"    Molecules: {len(molDf):,}")

    # InChIKey 碰撞 (InChIKey collision)
    glycoKeys = set(df["standard_inchi_key"].dropna().astype(str))
    chemblKeys = set(molDf["standard_inchi_key"].astype(str))
    hitKeys = glycoKeys & chemblKeys
    print(f"    InChIKey hits: {len(hitKeys):,} / {len(glycoKeys):,} "
          f"({len(hitKeys)/len(glycoKeys)*100:.1f}%)")

    if not hitKeys:
        print(f"    [WARN] No InChIKey hits!")
        return df

    # 提取命中的 ChEMBL IDs
    hitMols = molDf[molDf["standard_inchi_key"].isin(hitKeys)]
    hitChemblIds = set(hitMols["molecule_chembl_id"].astype(str))

    # 加载活性 — 仅命中 IDs (Load activities — only hit IDs)
    print(f"  [ChEMBL] Loading activities: {chemblActPath}")
    actChunks = []
    chunkSize = 500_000
    for chunk in pd.read_csv(chemblActPath, dtype=str, low_memory=False,
                             chunksize=chunkSize):
        filtered = chunk[chunk["molecule_chembl_id"].isin(hitChemblIds)]
        if len(filtered) > 0:
            actChunks.append(filtered)
        print(f"    ... scanned {chunkSize} rows, kept {len(filtered)}", end="\r")

    if not actChunks:
        print(f"\n    [WARN] No activity data for hit molecules")
        return df

    actDf = pd.concat(actChunks, ignore_index=True)
    print(f"\n    Activity rows for hits: {len(actDf):,}")

    # 合并: InChIKey → chembl_id → activities
    keyToChembl = hitMols.set_index("standard_inchi_key")["molecule_chembl_id"].to_dict()

    # 为每个 GlycoNP 化合物聚合其 ChEMBL 活性
    targetCol = "target_pref_name" if "target_pref_name" in actDf.columns else "target_name"
    typeCol = "standard_type" if "standard_type" in actDf.columns else "type"
    pchemblCol = "pchembl_value" if "pchembl_value" in actDf.columns else None

    df["ChEMBL_Targets"] = ""
    df["ChEMBL_Activity_Types"] = ""
    df["ChEMBL_Best_pChEMBL"] = np.nan

    filled = 0
    for idx in df.index:
        ik = str(df.at[idx, "standard_inchi_key"])
        if ik not in keyToChembl:
            continue

        chemblId = keyToChembl[ik]
        compAct = actDf[actDf["molecule_chembl_id"] == chemblId]

        if len(compAct) == 0:
            continue

        # 聚合靶点 (Aggregate targets)
        targets = compAct[targetCol].dropna().unique()
        types = compAct[typeCol].dropna().unique()

        df.at[idx, "ChEMBL_Targets"] = " | ".join(targets[:5])
        df.at[idx, "ChEMBL_Activity_Types"] = ", ".join(types[:5])

        if pchemblCol and pchemblCol in compAct.columns:
            pvals = pd.to_numeric(compAct[pchemblCol], errors="coerce")
            if pvals.notna().any():
                df.at[idx, "ChEMBL_Best_pChEMBL"] = pvals.max()

        filled += 1

    print(f"  [ChEMBL] Filled {filled:,} compounds with activity data")
    return df


# =====================================================================
# Engine 2: OpenAlex 合规文献概念提取 (OpenAlex Polite Pool Mining)
# =====================================================================

def mineOpenAlexConcepts(
    df: pd.DataFrame,
    pilotSize: int = 500,
    mailto: str = "glyconp_pipeline@example.com",
    minScore: float = 0.5,
    targetClasses: Optional[List[str]] = None,
) -> pd.DataFrame:
    """
    通过 OpenAlex Polite Pool API 提取文献概念关键词。
    Extract literature concepts via OpenAlex Polite Pool API.

    合规性措施 (Compliance Measures):
      1. mailto 参数进入 Polite Pool (合法高速通道)
      2. time.sleep(0.1) 每次请求间隔 100ms
      3. 仅提取 concepts 元数据, 不下载全文
      4. 先导测试仅 500 样本, 不会对服务器造成压力

    API 格式 (API Format):
      GET https://api.openalex.org/works/https://doi.org/{doi}?mailto={email}
      Response JSON → .concepts[] → [{display_name, score}]
    """
    import requests

    if targetClasses is None:
        targetClasses = ["Steroids", "Flavonoids", "Glycopeptide", "Glycolipids"]

    print(f"\n  [OpenAlex] Polite Pool Concept Mining")
    print(f"    mailto: {mailto}")
    print(f"    Target classes: {targetClasses}")
    print(f"    Min concept score: {minScore}")

    # 筛选高价值化合物: 有 DOI + 属于目标大类
    hasDoi = df["dois"].notna() & (~df["dois"].astype(str).isin(["", "nan", "None"]))

    # 清理 Superclass
    def cleanClass(val):
        s = str(val).strip()
        if "(Tanimoto=" in s:
            s = s[:s.index("(Tanimoto=")].strip()
        return s

    df["_class_clean"] = df["Superclass"].apply(cleanClass)
    inTarget = df["_class_clean"].isin(targetClasses)

    candidates = df[hasDoi & inTarget].copy()
    print(f"    Candidates (has DOI + target class): {len(candidates):,}")

    if len(candidates) == 0:
        print(f"    [WARN] No candidates found, expanding to all classes")
        candidates = df[hasDoi].copy()
        print(f"    Expanded candidates: {len(candidates):,}")

    # 随机抽样 (Random sampling)
    nSample = min(pilotSize, len(candidates))
    pilotDf = candidates.sample(n=nSample, random_state=42)
    print(f"    Pilot sample: {nSample}")

    # API 调用 (API calls)
    BASE_URL = "https://api.openalex.org/works/https://doi.org/"

    # 准备结果列
    df["OpenAlex_Concepts"] = ""

    success = 0
    failed = 0
    allConcepts: List[str] = []

    for i, (idx, row) in enumerate(pilotDf.iterrows()):
        doiRaw = str(row["dois"]).strip()
        # 取第一个 DOI (如果有多个, 用 | 分隔)
        doi = doiRaw.split("|")[0].strip()
        # 清理 DOI 格式
        doi = doi.replace("https://doi.org/", "").replace("http://doi.org/", "")
        doi = doi.strip("/").strip()

        if not doi or doi == "nan":
            failed += 1
            continue

        try:
            url = f"{BASE_URL}{doi}"
            # 合规 Polite Pool: mailto 参数 + User-Agent
            params = {"mailto": mailto}
            headers = {
                "User-Agent": f"GlycoNP-Pipeline/1.0 (mailto:{mailto})",
                "Accept": "application/json",
            }

            resp = requests.get(url, params=params, headers=headers, timeout=10)

            if resp.status_code == 200:
                data = resp.json()

                # 提取概念 (Extract concepts with score > threshold)
                concepts = data.get("concepts", [])
                # OpenAlex v2 可能用 topics 替代 concepts
                if not concepts:
                    topics = data.get("topics", [])
                    concepts = [{"display_name": t.get("display_name", ""),
                                 "score": t.get("score", 0)} for t in topics]

                highConcepts = [
                    c["display_name"] for c in concepts
                    if c.get("score", 0) >= minScore and c.get("display_name")
                ]

                if highConcepts:
                    df.at[idx, "OpenAlex_Concepts"] = " | ".join(highConcepts)
                    allConcepts.extend(highConcepts)
                    success += 1
                else:
                    success += 1  # API 成功但无高分概念

            elif resp.status_code == 404:
                failed += 1
            else:
                failed += 1

        except Exception as e:
            failed += 1

        # 进度显示 (Progress display)
        if (i + 1) % 50 == 0:
            print(f"    Progress: {i+1}/{nSample} "
                  f"(success={success}, failed={failed})")

        # 合规间隔 (Polite interval — 100ms)
        time.sleep(0.1)

    print(f"\n  [OpenAlex] Results:")
    print(f"    Success: {success}")
    print(f"    Failed/404: {failed}")
    print(f"    Total concepts extracted: {len(allConcepts)}")

    df.drop(columns=["_class_clean"], inplace=True, errors="ignore")
    return df, allConcepts


# =====================================================================
# 3. 糖链-活性映射表 (Sugar–Activity Mapping)
# =====================================================================

def buildSugarActivityMapping(
    df: pd.DataFrame,
    allConcepts: List[str],
    reportDir: str,
):
    """
    构建糖链修饰与生物活性的交叉映射。
    Build cross-mapping between sugar modifications and bioactivity concepts.
    """
    # Top 概念词频 (Top concept frequencies)
    conceptCounter = Counter(allConcepts)

    mdLines = ["# GlycoNP Bioactivity & Literature Mining Report",
               "# 糖缀合物生物活性与文献概念挖掘报告\n"]

    # ---- Part A: Top concepts ----
    mdLines.append("## Top 20 OpenAlex Concepts (score > 0.5)\n")
    mdLines.append("| Rank | Concept | Count |")
    mdLines.append("|:----:|:--------|------:|")

    print(f"\n  [Mapping] Top 20 concepts:")
    for i, (concept, count) in enumerate(conceptCounter.most_common(20), 1):
        print(f"    {i:2d}. {concept:40s} {count}")
        mdLines.append(f"| {i} | {concept} | {count} |")

    # ---- Part B: ChEMBL target summary ----
    hasTarget = df["ChEMBL_Targets"].notna() & (df["ChEMBL_Targets"] != "")
    targetRows = df[hasTarget]

    allTargets = []
    for targets in targetRows["ChEMBL_Targets"]:
        for t in str(targets).split(" | "):
            t = t.strip()
            if t and t != "nan":
                allTargets.append(t)

    targetCounter = Counter(allTargets)
    mdLines.append("\n## Top 20 ChEMBL Targets\n")
    mdLines.append("| Rank | Target | Compounds |")
    mdLines.append("|:----:|:-------|----------:|")

    print(f"\n  [Mapping] Top 20 ChEMBL targets:")
    for i, (target, count) in enumerate(targetCounter.most_common(20), 1):
        print(f"    {i:2d}. {target:45s} {count}")
        mdLines.append(f"| {i} | {target} | {count} |")

    # ---- Part C: Modification × Activity cross-mapping ----
    mdLines.append("\n## Sugar Modification → Activity Enrichment\n")
    mdLines.append(
        "> For each major modification, what are the most common "
        "activity concepts?\n")

    targetMods = ["O-Me", "O-Ac", "NAc", "COOH", "NH2", "Sulfate", "Phosphate"]

    hasConceptOrTarget = (
        (df["OpenAlex_Concepts"].fillna("") != "")
        | (df["ChEMBL_Targets"].fillna("") != "")
    )
    activeDf = df[hasConceptOrTarget].copy()

    if len(activeDf) > 0:
        mdLines.append("| Modification | N molecules | Top Concepts/Targets |")
        mdLines.append("|:-------------|:-----------:|:---------------------|")

        print(f"\n  [Mapping] Modification → Activity:")
        for mod in targetMods:
            modMask = activeDf["Glycan_Modifications"].fillna("").str.contains(
                mod, na=False)
            modDf = activeDf[modMask]
            nMol = len(modDf)

            if nMol == 0:
                continue

            # 收集该修饰的所有概念+靶点
            modConcepts = []
            for _, row in modDf.iterrows():
                concepts = str(row.get("OpenAlex_Concepts", ""))
                targets = str(row.get("ChEMBL_Targets", ""))
                for c in concepts.split(" | "):
                    if c.strip() and c.strip() != "nan":
                        modConcepts.append(c.strip())
                for t in targets.split(" | "):
                    if t.strip() and t.strip() != "nan":
                        modConcepts.append(t.strip())

            modCounter = Counter(modConcepts)
            topTerms = ", ".join([f"**{t}**({c})"
                                  for t, c in modCounter.most_common(3)])
            print(f"    {mod:12s} n={nMol:>4} → {topTerms}")
            mdLines.append(f"| `{mod}` | {nMol} | {topTerms} |")

    mdLines.append(
        "\n---\n> Generated by `mine_bioactivity_and_literature.py` | "
        "GlycoNP Pipeline Endgame")

    mdPath = os.path.join(reportDir, "Bioactivity_Literature_Report.md")
    with open(mdPath, "w", encoding="utf-8") as f:
        f.write("\n".join(mdLines))
    print(f"\n  Report saved: {mdPath}")

    return mdPath


# =====================================================================
# Main
# =====================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Bioactivity & Literature Mining")
    parser.add_argument("--input", type=str, default=None)
    parser.add_argument("--pilot", type=int, default=500,
                        help="OpenAlex pilot sample size")
    parser.add_argument("--email", type=str,
                        default="glyconp_pipeline@example.com",
                        help="Email for OpenAlex Polite Pool")
    args = parser.parse_args()

    baseDir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    reportDir = os.path.join(baseDir, "reports")
    dataDir = os.path.join(baseDir, "data")

    if args.input:
        inputPath = args.input
    else:
        inputPath = os.path.join(reportDir, "GlycoNP_Pipeline_Full_Cleaned.csv")

    print("=" * 70)
    print("  GlycoNP Bioactivity & Literature Mining — Endgame")
    print("  糖缀合物生物活性与文献概念挖掘")
    print("=" * 70)
    print(f"  Input: {inputPath}")

    t0 = time.time()
    df = pd.read_csv(inputPath, low_memory=False, dtype=str,
                     encoding="utf-8-sig")
    print(f"  Loaded: {len(df):,} rows ({time.time()-t0:.1f}s)")

    # ---- Engine 1: ChEMBL ----
    print(f"\n{'='*70}")
    print(f"  Engine 1: ChEMBL Local Static Merge")
    print(f"{'='*70}")

    chemblMolPath = os.path.join(dataDir, "chembl", "chembl_molecules.csv")
    chemblActPath = os.path.join(dataDir, "chembl", "chembl_activities.csv")

    if os.path.exists(chemblMolPath) and os.path.exists(chemblActPath):
        df = mergeChembl(df, chemblMolPath, chemblActPath)
    else:
        print(f"  [SKIP] ChEMBL CSVs not found:")
        print(f"    Molecules: {chemblMolPath}")
        print(f"    Activities: {chemblActPath}")

    # ---- Engine 2: OpenAlex ----
    print(f"\n{'='*70}")
    print(f"  Engine 2: OpenAlex Polite Pool Concept Mining")
    print(f"{'='*70}")

    df, allConcepts = mineOpenAlexConcepts(
        df,
        pilotSize=args.pilot,
        mailto=args.email,
    )

    # ---- Mapping ----
    print(f"\n{'='*70}")
    print(f"  Sugar–Activity Mapping")
    print(f"{'='*70}")

    mdPath = buildSugarActivityMapping(df, allConcepts, reportDir)

    # ---- Save ----
    outPath = os.path.join(reportDir, "GlycoNP_Bioactivity_Mined.csv")
    df.to_csv(outPath, index=False, encoding="utf-8-sig")

    elapsed = time.time() - t0
    print(f"\n{'='*70}")
    print(f"  Endgame Complete!")
    print(f"  Time: {elapsed:.0f}s")
    print(f"  Output: {outPath}")
    print(f"  Report: {mdPath}")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
