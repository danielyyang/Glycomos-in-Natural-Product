"""
GlycoNP × ChEMBL 生物活性集成引擎
GlycoNP × ChEMBL Bioactivity Integration Engine

将管线输出的糖缀合物通过 InChIKey 与 ChEMBL 数据库进行碰撞,
提取已知的生物靶点和活性数据, 揭示"糖链特征 → 生物靶点"的关联。

工作流 (Workflow):
  1. 加载 GlycoNP 全量管线结果 (InChIKey, 糖链序列, 修饰, 大类...)
  2. 加载 ChEMBL 分子表 (InChIKey → molecule_chembl_id)
  3. 加载 ChEMBL 活性表 (chembl_id → 靶点, IC50, Ki, pChEMBL)
  4. 三表 Merge → 生成 GlycoNP_Final_Discovery_DB.csv
  5. "活性-糖链"相关性分析报告

使用方法 (Usage):
  # 真实模式 (Real mode) — 需要 ChEMBL 文件
  python scripts/integrate_bioactivity.py

  # 演示模式 (Mock mode) — 模拟 500 个命中, 无需真实文件
  python scripts/integrate_bioactivity.py --mock

ChEMBL 文件准备 (Data Preparation):
  用户需将以下文件放入 data/chembl/ 目录:
  1. chembl_molecules.csv — 分子映射表
     必须包含列: standard_inchi_key, molecule_chembl_id
  2. chembl_activities.csv — 活性数据表
     必须包含列: molecule_chembl_id, target_pref_name,
                  standard_type, standard_value, standard_units,
                  pchembl_value
"""
import argparse
import os
import sys
import random
import hashlib
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))


# =====================================================================
# 1. 数据加载 (Data Loading)
# =====================================================================

def loadGlycoNpData(path: str) -> pd.DataFrame:
    """
    加载 GlycoNP 管线全量输出。
    Load GlycoNP pipeline full output.
    """
    print(f"  [GlycoNP] Loading: {path}")
    df = pd.read_csv(path, low_memory=False, dtype=str, encoding="utf-8-sig")
    print(f"    Rows: {len(df):,} | InChIKeys: {df['standard_inchi_key'].nunique():,}")
    return df


def loadChemblMolecules(path: str) -> pd.DataFrame:
    """
    加载 ChEMBL 分子映射表 (InChIKey → chembl_id)。
    Load ChEMBL molecule lookup table.

    只保留必要列以节省内存 (Keep only essential columns).
    ChEMBL 分子表通常有 200 万行, 需要高效处理。
    """
    print(f"  [ChEMBL] Loading molecules: {path}")
    # 只读取必要列 (Read only needed columns)
    useCols = ["standard_inchi_key", "molecule_chembl_id"]
    df = pd.read_csv(path, usecols=useCols, dtype=str, low_memory=False)
    df = df.dropna(subset=["standard_inchi_key"]).drop_duplicates()
    print(f"    Unique molecules: {len(df):,}")
    return df


def loadChemblActivities(
    path: str,
    chemblIds: Optional[set] = None,
) -> pd.DataFrame:
    """
    加载 ChEMBL 活性数据并过滤。
    Load ChEMBL activity data, optionally filtering by molecule IDs.

    设计意图: ChEMBL 活性表有数千万行。先做 InChIKey 碰撞筛出命中分子的
    chembl_id 列表, 再只保留这些分子的活性数据, 大幅减少内存占用。
    """
    print(f"  [ChEMBL] Loading activities: {path}")
    useCols = [
        "molecule_chembl_id",
        "target_pref_name",
        "target_organism",
        "standard_type",
        "standard_value",
        "standard_units",
        "pchembl_value",
        "assay_type",
    ]
    # 分块读取以控制内存 (Chunked reading for memory control)
    chunks = []
    chunkSize = 500_000
    totalRead = 0
    totalKept = 0

    for chunk in pd.read_csv(
        path, usecols=useCols, dtype=str,
        low_memory=False, chunksize=chunkSize,
    ):
        totalRead += len(chunk)
        if chemblIds is not None:
            chunk = chunk[chunk["molecule_chembl_id"].isin(chemblIds)]
        totalKept += len(chunk)
        if not chunk.empty:
            chunks.append(chunk)

    if chunks:
        df = pd.concat(chunks, ignore_index=True)
    else:
        df = pd.DataFrame(columns=useCols)

    print(f"    Read {totalRead:,} rows -> kept {totalKept:,} for hit molecules")
    return df


# =====================================================================
# 2. 数据碰撞与合并 (Merge / Data Collision)
# =====================================================================

def mergeGlycoNpWithChembl(
    glycoDf: pd.DataFrame,
    chemblMolDf: pd.DataFrame,
    chemblActDf: pd.DataFrame,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    三表合并: GlycoNP ←[InChIKey]→ ChEMBL_mol ←[chembl_id]→ ChEMBL_act.
    Three-table merge via InChIKey and molecule_chembl_id.

    Returns:
        (discoveryDf, hitSummaryDf):
          discoveryDf — 完整的发现数据库 (每行 = 一个 化合物×活性 记录)
          hitSummaryDf — 命中化合物汇总 (每行 = 一个唯一化合物)
    """
    print("\n  [Merge] Step 1: InChIKey collision...")
    # Left join: 保留所有 GlycoNP 化合物 (Keep all GlycoNP compounds)
    merged = glycoDf.merge(
        chemblMolDf, on="standard_inchi_key", how="left",
    )
    hitCount = merged["molecule_chembl_id"].notna().sum()
    print(f"    GlycoNP compounds with ChEMBL match: {hitCount:,} / {len(glycoDf):,}"
          f" ({hitCount/len(glycoDf)*100:.1f}%)")

    # 提取命中的 chembl_id (Extract hit chembl_ids)
    hitDf = merged[merged["molecule_chembl_id"].notna()].copy()

    if chemblActDf.empty:
        print("    [WARNING] No activity data available")
        return merged, hitDf

    print("  [Merge] Step 2: Activity data join...")
    # 将活性数据接入命中化合物 (Join activities to hit compounds)
    discoveryDf = hitDf.merge(
        chemblActDf, on="molecule_chembl_id", how="left",
    )
    withActivity = discoveryDf["target_pref_name"].notna().sum()
    print(f"    Compound-activity pairs: {withActivity:,}")
    print(f"    Unique targets: {discoveryDf['target_pref_name'].dropna().nunique():,}")

    return discoveryDf, hitDf


# =====================================================================
# 3. 活性-糖链相关性分析 (Activity-Glycan Correlation Analysis)
# =====================================================================

def analyzeActivityGlycanCorrelation(
    discoveryDf: pd.DataFrame,
    outputDir: str,
) -> str:
    """
    核心分析: 糖链特征 → 生物靶点 的关联模式。
    Core analysis: Sugar features → Biological target associations.

    分析维度 (Analysis Dimensions):
      A. 修饰-靶点关联 (Modification → Target)
      B. 糖链序列-靶点关联 (Sugar Sequence → Target)
      C. 高活性化合物的糖链模式 (High-potency glycan patterns)
      D. 大类-活性分布 (Superclass → Activity distribution)
    """
    lines = []
    lines.append("# GlycoNP Discovery Report: Sugar-Activity Correlation")
    lines.append("# GlycoNP 发现报告: 糖链-活性关联分析\n")
    lines.append(f"> Total compound-activity pairs: **{len(discoveryDf):,}**")
    lines.append(f"> Unique targets: "
                 f"**{discoveryDf['target_pref_name'].dropna().nunique():,}**\n")

    # ------------------------------------------------------------------
    # A. 修饰-靶点关联 (Modification → Target Association)
    # ------------------------------------------------------------------
    lines.append("## A. Modification → Target Association")
    lines.append("## A. 修饰基团-靶点关联\n")

    modCols = ["O-Me", "O-Ac", "NAc", "Sulfate", "Phosphate", "COOH", "NH2"]
    modTargetRows = []

    for _, row in discoveryDf.iterrows():
        mods = str(row.get("Glycan_Modifications", ""))
        target = str(row.get("target_pref_name", ""))
        if target == "nan" or not target:
            continue
        for mod in modCols:
            if mod in mods:
                modTargetRows.append({"modification": mod, "target": target})

    if modTargetRows:
        modTargetDf = pd.DataFrame(modTargetRows)
        lines.append("| Modification | Top Target | Compound Count |")
        lines.append("|:------------|:-----------|---------------:|")

        for mod in modCols:
            subset = modTargetDf[modTargetDf["modification"] == mod]
            if subset.empty:
                continue
            topTarget = subset["target"].value_counts().head(1)
            if not topTarget.empty:
                lines.append(
                    f"| **{mod}** | {topTarget.index[0][:50]} | "
                    f"{int(topTarget.values[0])} |"
                )

        lines.append("")

        # 详细: 每种修饰的 Top 3 靶点 (Detailed: Top 3 targets per mod)
        lines.append("### Detailed: Top 3 Targets per Modification\n")
        for mod in modCols:
            subset = modTargetDf[modTargetDf["modification"] == mod]
            if len(subset) < 3:
                continue
            lines.append(f"**{mod}** ({len(subset)} hit compound-target pairs):\n")
            lines.append("| Rank | Target | Count |")
            lines.append("|:----:|--------|------:|")
            for i, (tgt, cnt) in enumerate(
                subset["target"].value_counts().head(3).items(), 1
            ):
                lines.append(f"| {i} | {tgt[:60]} | {int(cnt)} |")
            lines.append("")

    # ------------------------------------------------------------------
    # B. 糖链序列-靶点关联 (Sugar Sequence → Target)
    # ------------------------------------------------------------------
    lines.append("## B. Sugar Sequence → Target Association")
    lines.append("## B. 糖链序列-靶点关联\n")

    validAct = discoveryDf[discoveryDf["target_pref_name"].notna()].copy()
    if not validAct.empty:
        # Top 10 糖链 × 靶点交叉 (Top sugar sequences × targets)
        topSugars = validAct["Sugar_Sequence"].value_counts().head(10).index.tolist()
        lines.append("| Sugar Sequence | #1 Target | #2 Target | Total Pairs |")
        lines.append("|:--------------|:----------|:----------|------------:|")

        for sugar in topSugars:
            sSubset = validAct[validAct["Sugar_Sequence"] == sugar]
            topTargets = sSubset["target_pref_name"].value_counts().head(2)
            t1 = topTargets.index[0][:30] if len(topTargets) > 0 else "-"
            t2 = topTargets.index[1][:30] if len(topTargets) > 1 else "-"
            lines.append(
                f"| `{sugar[:35]}` | {t1} | {t2} | {len(sSubset)} |"
            )
        lines.append("")

    # ------------------------------------------------------------------
    # C. 高活性化合物的糖链模式 (High-potency patterns)
    # ------------------------------------------------------------------
    lines.append("## C. High-Potency Compound Glycan Patterns (pChEMBL >= 7)")
    lines.append("## C. 高活性化合物的糖链模式 (pChEMBL >= 7, 即 IC50 < 100nM)\n")

    discoveryDf["_pchembl"] = pd.to_numeric(
        discoveryDf.get("pchembl_value", pd.Series(dtype="float64")),
        errors="coerce",
    )
    potent = discoveryDf[discoveryDf["_pchembl"] >= 7.0]

    if not potent.empty:
        lines.append(f"Compounds with pChEMBL >= 7: **{len(potent)}**\n")

        lines.append("| Sugar Sequence | Target | pChEMBL | Modifications |")
        lines.append("|:-------------|:-------|--------:|:-------------|")

        for _, row in potent.nlargest(15, "_pchembl").iterrows():
            sugar = str(row.get("Sugar_Sequence", ""))[:30]
            target = str(row.get("target_pref_name", ""))[:35]
            pval = row["_pchembl"]
            mods = str(row.get("Glycan_Modifications", ""))
            if mods == "nan":
                mods = "None"
            lines.append(f"| `{sugar}` | {target} | {pval:.1f} | {mods} |")
        lines.append("")
    else:
        lines.append("No compounds with pChEMBL >= 7 found.\n")

    # ------------------------------------------------------------------
    # D. 大类-活性分布 (Superclass → Activity)
    # ------------------------------------------------------------------
    lines.append("## D. Superclass Activity Distribution")
    lines.append("## D. 大类活性分布\n")

    if not validAct.empty:
        # 清理 Superclass (Clean superclass)
        validAct["_class"] = validAct["Superclass"].apply(
            lambda x: _cleanClass(str(x)))

        classAct = validAct.groupby("_class")["target_pref_name"].nunique()
        classAct = classAct.sort_values(ascending=False).head(10)

        lines.append("| Superclass | Unique Targets | Total Pairs |")
        lines.append("|:-----------|---------------:|------------:|")

        for cls, nTargets in classAct.items():
            nPairs = len(validAct[validAct["_class"] == cls])
            lines.append(f"| {cls} | {int(nTargets)} | {nPairs} |")
        lines.append("")

    # ------------------------------------------------------------------
    # E. 终极发现 (Key Discoveries)
    # ------------------------------------------------------------------
    lines.append("## E. Key Discoveries / 终极发现\n")

    # 找出修饰基团与靶点的最强关联
    if modTargetRows:
        modTargetDf = pd.DataFrame(modTargetRows)
        # 按修饰-靶点对排序
        topPairs = modTargetDf.groupby(
            ["modification", "target"]
        ).size().sort_values(ascending=False).head(5)

        lines.append("### Top 5 Modification-Target Pairs\n")
        lines.append("| Rank | Modification | Target | Count |")
        lines.append("|:----:|:------------|:-------|------:|")
        for i, ((mod, tgt), cnt) in enumerate(topPairs.items(), 1):
            lines.append(f"| {i} | **{mod}** | {tgt[:50]} | {int(cnt)} |")
        lines.append("")

    lines.append("---")
    lines.append("> Generated by `integrate_bioactivity.py` | GlycoNP Pipeline")

    # 保存报告 (Save report)
    reportPath = os.path.join(outputDir, "Discovery_Activity_Report.md")
    with open(reportPath, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))

    print(f"\n  Report saved: {reportPath}")
    return reportPath


def _cleanClass(val: str) -> str:
    """清理 Superclass 名称 (Clean Superclass)."""
    if not val or val == "nan":
        return "Unclassified"
    s = str(val).strip()
    if "(Tanimoto=" in s:
        s = s[:s.index("(Tanimoto=")].strip()
    if s.startswith("Glycolipid"):
        s = "Glycolipids"
    return s if s else "Unclassified"


# =====================================================================
# 4. Mock 数据生成 (Mock Data Generator)
# =====================================================================

def generateMockChemblData(
    glycoDf: pd.DataFrame,
    nHits: int = 500,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    生成模拟的 ChEMBL 碰撞数据, 用于测试分析逻辑。
    Generate simulated ChEMBL hit data for testing the analysis pipeline.

    设计原则: 模拟数据应尽可能反映真实分布。
    不同大类配置不同的靶点概率, 模拟真实的药理学偏好。
    注意: 所有数值为合成数据 [TEST DATA ONLY]。
    """
    print(f"\n  [MOCK] Generating simulated ChEMBL data ({nHits} hits)...")

    random.seed(42)
    np.random.seed(42)

    # 靶点定义: 按治疗领域分组 (Targets by therapeutic area)
    # [TEST DATA ONLY] — 靶点名称来自公开文献
    TARGET_POOLS = {
        "anti_inflammatory": [
            "Cyclooxygenase-2 (COX-2)",
            "5-lipoxygenase (5-LOX)",
            "Tumor necrosis factor alpha (TNF-alpha)",
            "Nuclear factor NF-kappa-B p65",
            "Interleukin-6 (IL-6)",
            "Phospholipase A2",
        ],
        "anti_tumor": [
            "Epidermal growth factor receptor (EGFR)",
            "Vascular endothelial growth factor receptor 2 (VEGFR2)",
            "Topoisomerase II alpha",
            "Tubulin beta chain",
            "Caspase-3",
            "Cyclin-dependent kinase 2 (CDK2)",
            "Poly [ADP-ribose] polymerase-1 (PARP-1)",
        ],
        "anti_microbial": [
            "DNA gyrase subunit B",
            "Dihydrofolate reductase (DHFR)",
            "Penicillin-binding protein 2a (PBP2a)",
            "Beta-lactamase",
        ],
        "metabolic": [
            "Alpha-glucosidase",
            "Dipeptidyl peptidase IV (DPP-4)",
            "Peroxisome proliferator-activated receptor gamma (PPAR-gamma)",
            "HMG-CoA reductase",
            "Aldose reductase",
        ],
        "neuro": [
            "Acetylcholinesterase (AChE)",
            "Monoamine oxidase B (MAO-B)",
            "GABA-A receptor",
            "Serotonin 5-HT2A receptor",
        ],
    }

    # 大类 → 靶点概率偏好 (Superclass → target area probability)
    # 反映真实药理学趋势: 黄酮类多抗炎, 生物碱多神经活性...
    CLASS_TARGET_BIAS = {
        "Flavonoids": {"anti_inflammatory": 0.4, "anti_tumor": 0.3,
                       "metabolic": 0.2, "neuro": 0.1},
        "Alkaloids": {"neuro": 0.4, "anti_tumor": 0.3,
                      "anti_microbial": 0.2, "anti_inflammatory": 0.1},
        "Phenylpropanoids": {"anti_inflammatory": 0.35, "anti_tumor": 0.25,
                             "metabolic": 0.25, "anti_microbial": 0.15},
        "Anthraquinones": {"anti_tumor": 0.4, "anti_microbial": 0.3,
                           "anti_inflammatory": 0.2, "metabolic": 0.1},
        "Coumarins": {"anti_inflammatory": 0.3, "anti_tumor": 0.3,
                      "metabolic": 0.2, "neuro": 0.2},
        "Steroids": {"anti_inflammatory": 0.45, "metabolic": 0.3,
                     "anti_tumor": 0.15, "neuro": 0.1},
        "Glycolipids": {"anti_microbial": 0.4, "anti_tumor": 0.3,
                        "anti_inflammatory": 0.2, "metabolic": 0.1},
        "Triterpenoids": {"anti_inflammatory": 0.35, "anti_tumor": 0.35,
                          "metabolic": 0.2, "anti_microbial": 0.1},
    }

    DEFAULT_BIAS = {"anti_inflammatory": 0.25, "anti_tumor": 0.25,
                    "metabolic": 0.2, "anti_microbial": 0.15, "neuro": 0.15}

    # 采样命中化合物 (Sample hit compounds)
    # 偏好有修饰的化合物 (Bias towards modified compounds)
    hasMods = glycoDf[glycoDf["Glycan_Modifications"].notna()
                      & (glycoDf["Glycan_Modifications"] != "")]
    noMods = glycoDf[~glycoDf.index.isin(hasMods.index)]

    nFromMods = min(int(nHits * 0.6), len(hasMods))
    nFromNoMods = min(nHits - nFromMods, len(noMods))

    hitRows = pd.concat([
        hasMods.sample(n=nFromMods, random_state=42) if nFromMods > 0 else pd.DataFrame(),
        noMods.sample(n=nFromNoMods, random_state=42) if nFromNoMods > 0 else pd.DataFrame(),
    ], ignore_index=True)

    # 为每个命中化合物生成 chembl_id 和活性数据
    mockMolRows = []
    mockActRows = []

    for idx, row in hitRows.iterrows():
        inchiKey = str(row.get("standard_inchi_key", ""))
        if not inchiKey or inchiKey == "nan":
            continue

        # 生成假 chembl_id (Generate mock chembl_id)
        hashVal = hashlib.md5(inchiKey.encode()).hexdigest()[:6]
        chemblId = f"CHEMBL{int(hashVal, 16) % 9999999}"

        mockMolRows.append({
            "standard_inchi_key": inchiKey,
            "molecule_chembl_id": chemblId,
        })

        # 确定靶点偏好 (Determine target bias by superclass)
        superclass = _cleanClass(str(row.get("Superclass", "")))
        bias = CLASS_TARGET_BIAS.get(superclass, DEFAULT_BIAS)

        # 每个化合物 1-3 个活性数据 (1-3 activities per compound)
        nActivities = random.choices([1, 2, 3], weights=[0.5, 0.35, 0.15])[0]

        for _ in range(nActivities):
            # 选择靶点领域 (Select target area)
            area = random.choices(
                list(bias.keys()), weights=list(bias.values())
            )[0]
            target = random.choice(TARGET_POOLS[area])

            # 生成活性值 (Generate activity values)
            # pChEMBL 正态分布: 均值 5.5 (IC50~3μM), 标准差 1.2
            pchembl = max(3.0, min(10.0, np.random.normal(5.5, 1.2)))
            ic50Nm = 10 ** (9 - pchembl)  # 换算 nM

            # 修饰基团会略微提升活性 (Modifications slightly boost activity)
            mods = str(row.get("Glycan_Modifications", ""))
            if "NAc" in mods or "Sulfate" in mods:
                pchembl += random.uniform(0.2, 0.8)
                ic50Nm = 10 ** (9 - pchembl)

            stdType = random.choice(["IC50", "Ki", "EC50"])

            mockActRows.append({
                "molecule_chembl_id": chemblId,
                "target_pref_name": target,
                "target_organism": "Homo sapiens",
                "standard_type": stdType,
                "standard_value": f"{ic50Nm:.1f}",
                "standard_units": "nM",
                "pchembl_value": f"{pchembl:.2f}",
                "assay_type": "B",
            })

    mockMolDf = pd.DataFrame(mockMolRows)
    mockActDf = pd.DataFrame(mockActRows)

    print(f"    Mock molecules: {len(mockMolDf)}")
    print(f"    Mock activities: {len(mockActDf)}")
    print(f"    Mock targets: {mockActDf['target_pref_name'].nunique()}")

    return mockMolDf, mockActDf


# =====================================================================
# 5. 主流程 (Main)
# =====================================================================

def main():
    parser = argparse.ArgumentParser(
        description="GlycoNP x ChEMBL Bioactivity Integration")
    parser.add_argument("--input", type=str, default=None,
                        help="GlycoNP pipeline CSV")
    parser.add_argument("--chembl-mol", type=str, default=None,
                        help="ChEMBL molecules CSV")
    parser.add_argument("--chembl-act", type=str, default=None,
                        help="ChEMBL activities CSV")
    parser.add_argument("--mock", action="store_true",
                        help="Use mock ChEMBL data for testing")
    parser.add_argument("--mock-hits", type=int, default=500,
                        help="Number of mock hits to simulate")
    parser.add_argument("--output", type=str, default=None,
                        help="Output discovery DB CSV")
    args = parser.parse_args()

    baseDir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    reportDir = os.path.join(baseDir, "reports")
    chemblDir = os.path.join(baseDir, "data", "chembl")

    # 确定输入路径 (Resolve input paths)
    if args.input:
        glycoPath = args.input
    else:
        glycoPath = os.path.join(reportDir, "GlycoNP_Pipeline_Full.csv")
        if not os.path.exists(glycoPath):
            glycoPath = os.path.join(reportDir, "GlycoNP_Pipeline_1000.csv")

    outputPath = args.output or os.path.join(
        reportDir, "GlycoNP_Final_Discovery_DB.csv")

    # Banner
    print("=" * 70)
    print("  GlycoNP x ChEMBL Bioactivity Integration Engine")
    print("  糖缀合物 × ChEMBL 生物活性集成引擎")
    print("=" * 70)

    # 加载 GlycoNP 数据 (Load GlycoNP data)
    glycoDf = loadGlycoNpData(glycoPath)

    if args.mock:
        # ============ MOCK MODE ============
        print("\n  *** MOCK MODE — Using simulated ChEMBL data ***")
        print("  *** [TEST DATA ONLY] ***\n")
        mockMolDf, mockActDf = generateMockChemblData(glycoDf, nHits=args.mock_hits)

        discoveryDf, hitDf = mergeGlycoNpWithChembl(glycoDf, mockMolDf, mockActDf)

    else:
        # ============ REAL MODE ============
        chemblMolPath = args.chembl_mol or os.path.join(
            chemblDir, "chembl_molecules.csv")
        chemblActPath = args.chembl_act or os.path.join(
            chemblDir, "chembl_activities.csv")

        if not os.path.exists(chemblMolPath):
            print(f"\n  [ERROR] ChEMBL molecule file not found: {chemblMolPath}")
            print(f"  Please download from https://www.ebi.ac.uk/chembl/")
            print(f"  and place in: {chemblDir}")
            print(f"\n  Or use --mock flag for testing.\n")
            return

        if not os.path.exists(chemblActPath):
            print(f"\n  [ERROR] ChEMBL activity file not found: {chemblActPath}")
            print(f"  Please download from https://www.ebi.ac.uk/chembl/")
            print(f"  and place in: {chemblDir}")
            print(f"\n  Or use --mock flag for testing.\n")
            return

        chemblMolDf = loadChemblMolecules(chemblMolPath)
        # 先碰撞获取命中 chembl_id (Collision first to get hit IDs)
        hitIds = set(
            glycoDf.merge(chemblMolDf, on="standard_inchi_key")
            ["molecule_chembl_id"].unique()
        )
        print(f"    InChIKey collision hits: {len(hitIds):,}")

        chemblActDf = loadChemblActivities(chemblActPath, chemblIds=hitIds)
        discoveryDf, hitDf = mergeGlycoNpWithChembl(
            glycoDf, chemblMolDf, chemblActDf)

    # 保存终极数据库 (Save discovery database)
    print(f"\n  [Save] Writing discovery DB: {outputPath}")
    discoveryDf.to_csv(outputPath, index=False, encoding="utf-8-sig")
    print(f"    Rows: {len(discoveryDf):,}")

    # 相关性分析报告 (Correlation analysis report)
    analyzeActivityGlycanCorrelation(discoveryDf, reportDir)

    # Summary
    print(f"\n{'='*70}")
    print(f"  Integration Complete!")
    print(f"{'='*70}")
    print(f"  Discovery DB: {outputPath}")
    print(f"  Rows: {len(discoveryDf):,}")
    targetCol = "target_pref_name"
    if targetCol in discoveryDf.columns:
        nTargets = discoveryDf[targetCol].dropna().nunique()
        print(f"  Unique targets: {nTargets}")
    print(f"{'='*70}\n")


if __name__ == "__main__":
    main()
