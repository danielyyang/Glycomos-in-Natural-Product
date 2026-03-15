"""
GlycoNP 糖链频率与修饰分析
GlycoNP Glycan Frequency & Modification Analysis

读取全量管线输出 CSV, 执行三层统计分析并生成 Markdown 报告:
  1. 按 Superclass 分组, Top 5 Sugar_Sequence 频率
  2. Glycolipid / Glycopeptide 特别修饰统计
  3. 物种来源 × 糖链交叉分析

使用方法 (Usage):
  python scripts/analyze_glycan_frequencies.py [--input PATH] [--output PATH]
"""
import argparse
import os
import sys
from collections import Counter
from typing import Dict, List, Tuple

import pandas as pd

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))


# =====================================================================
# 1. 糖链频率分析 (Sugar Sequence Frequency by Superclass)
# =====================================================================

def analyzeSugarByClass(df: pd.DataFrame, topN: int = 5) -> Dict[str, List[Tuple[str, int]]]:
    """
    按 Superclass 分组, 统计每个大类下 Sugar_Sequence 的频率 Top N.
    Group by Superclass and count top N Sugar_Sequence per class.

    Returns:
        {superclass: [(sequence, count), ...]}
    """
    results = {}
    validDf = df[df["Sugar_Sequence"].notna() & (df["Sugar_Sequence"] != "") & (df["Sugar_Sequence"] != "nan")]

    for cls, group in validDf.groupby("Superclass"):
        if not cls or str(cls) == "nan":
            continue
        counter = Counter(group["Sugar_Sequence"].astype(str))
        results[str(cls)] = counter.most_common(topN)

    return results


# =====================================================================
# 2. 特定大类修饰统计 (Modification Stats for Special Classes)
# =====================================================================

def analyzeModsByClass(df: pd.DataFrame, targetClasses: List[str]) -> Dict[str, Counter]:
    """
    针对指定大类, 统计其 Glycan_Modifications 中各修饰基团的出现频率.
    For target classes, count frequency of each modification type.

    Returns:
        {class_name: Counter({mod_name: count, ...})}
    """
    results = {}
    validDf = df[
        df["Glycan_Modifications"].notna()
        & (df["Glycan_Modifications"] != "")
        & (df["Glycan_Modifications"] != "nan")
    ]

    for cls in targetClasses:
        classDf = validDf[validDf["Superclass"] == cls]
        if classDf.empty:
            # 尝试模糊匹配 (Fuzzy match for Tanimoto-appended names)
            classDf = validDf[validDf["Superclass"].str.contains(cls, case=False, na=False)]

        modCounter = Counter()
        for modStr in classDf["Glycan_Modifications"]:
            # 格式: "O-Ac: 2, NAc: 1"
            for part in str(modStr).split(","):
                part = part.strip()
                if ":" in part:
                    modName = part.split(":")[0].strip()
                    try:
                        modCount = int(part.split(":")[1].strip())
                    except ValueError:
                        modCount = 1
                    modCounter[modName] += modCount

        results[cls] = modCounter

    return results


# =====================================================================
# 3. 物种来源交叉分析 (Organism × Glycan Cross-Analysis)
# =====================================================================

def analyzeOrganismGlycan(df: pd.DataFrame, topN: int = 10) -> List[Tuple[str, str, int]]:
    """
    统计最常见的 (Organism, Sugar_Sequence) 组合.
    Count most common (Organism, Sugar_Sequence) pairs.

    Returns:
        [(organism, sugar_sequence, count), ...]
    """
    validDf = df[
        df["organisms"].notna() & (df["organisms"] != "nan")
        & df["Sugar_Sequence"].notna() & (df["Sugar_Sequence"] != "nan")
    ].copy()

    if validDf.empty:
        return []

    # 拆分多物种行 (Explode multi-organism rows: "Org1|Org2" → 2 rows)
    validDf["organisms"] = validDf["organisms"].astype(str)
    pairs = Counter()
    for _, row in validDf.iterrows():
        orgs = str(row["organisms"]).split("|")
        seq = str(row["Sugar_Sequence"])
        for org in orgs:
            org = org.strip()
            if org and org != "nan":
                pairs[(org, seq)] += 1

    return [(org, seq, count) for (org, seq), count in pairs.most_common(topN)]


# =====================================================================
# 4. Markdown 报告生成 (Report Generation)
# =====================================================================

def generateReport(
    df: pd.DataFrame,
    outputPath: str,
    topNSeq: int = 5,
    topNOrg: int = 15,
):
    """
    生成 Glycomics 汇总报告 (Markdown).
    Generate Glycomics Summary Report.
    """
    lines = []

    # Header
    lines.append("# GlycoNP 糖链组学汇总报告")
    lines.append("# GlycoNP Glycomics Summary Report\n")
    lines.append(f"> 数据来源: COCONUT 糖苷子集 | 总化合物: **{len(df):,}**\n")

    # Overall stats
    lines.append("## 一、全局统计 (Global Statistics)\n")

    nucCount = df["Has_Nucleotide"].sum() if "Has_Nucleotide" in df.columns else 0
    pepCount = df["Has_Peptide"].sum() if "Has_Peptide" in df.columns else 0
    modCount = sum(1 for v in df.get("Glycan_Modifications", []) if v and str(v) not in ("", "nan"))
    seqCount = sum(1 for v in df.get("Sugar_Sequence", []) if v and str(v) not in ("", "nan"))

    lines.append("| 指标 | 数值 |")
    lines.append("|------|------|")
    lines.append(f"| 总化合物数 | {len(df):,} |")
    lines.append(f"| 含糖序列 | {seqCount:,} ({seqCount/len(df)*100:.1f}%) |")
    lines.append(f"| 含修饰基团 | {modCount:,} ({modCount/len(df)*100:.1f}%) |")
    lines.append(f"| 核苷酸糖 | {nucCount:,} |")
    lines.append(f"| 糖肽/氨基酸糖苷 | {pepCount:,} |")

    # Superclass distribution
    lines.append("\n## 二、大类分布 (Superclass Distribution)\n")
    classDist = df["Superclass"].value_counts()
    lines.append("| Superclass | Count | % |")
    lines.append("|-----------|------:|----:|")
    for cls, count in classDist.head(15).items():
        pct = count / len(df) * 100
        lines.append(f"| {cls} | {count:,} | {pct:.1f}% |")

    # Sugar sequence by class
    lines.append("\n## 三、各大类 Top 5 糖链序列 (Sugar Sequences by Superclass)\n")
    sugarByClass = analyzeSugarByClass(df, topN=topNSeq)

    # 按化合物总数从多到少排列 (Sort by class size descending)
    sortedClasses = sorted(sugarByClass.keys(),
                           key=lambda c: classDist.get(c, 0), reverse=True)

    for cls in sortedClasses[:20]:  # Top 20 classes
        seqs = sugarByClass[cls]
        if not seqs:
            continue
        classSize = classDist.get(cls, 0)
        lines.append(f"### {cls} (n={classSize:,})\n")
        lines.append("| Rank | Sugar Sequence | Count | % of Class |")
        lines.append("|:----:|----------------|------:|-----------:|")
        for rank, (seq, count) in enumerate(seqs, 1):
            pct = count / classSize * 100 if classSize > 0 else 0
            lines.append(f"| {rank} | `{seq}` | {count:,} | {pct:.1f}% |")
        lines.append("")

    # Glycolipid / Glycopeptide modification analysis
    lines.append("\n## 四、特殊大类修饰统计 (Modifications in Special Classes)\n")
    targetClasses = ["Glycolipid", "Glycopeptide", "Nucleotide Sugar", "Flavonoids", "Steroids"]
    modsByClass = analyzeModsByClass(df, targetClasses)

    for cls, modCounter in modsByClass.items():
        if not modCounter:
            lines.append(f"### {cls}\n")
            lines.append("> 无修饰基团命中\n")
            continue
        total = sum(modCounter.values())
        lines.append(f"### {cls} (总修饰命中: {total})\n")
        lines.append("| Modification | Count | % |")
        lines.append("|-------------|------:|----:|")
        for modName, count in modCounter.most_common():
            pct = count / total * 100
            lines.append(f"| {modName} | {count:,} | {pct:.1f}% |")
        lines.append("")

    # Organism × Glycan cross-analysis
    lines.append("\n## 五、物种-糖链交叉分析 (Organism × Sugar Sequence)\n")
    orgGlycan = analyzeOrganismGlycan(df, topN=topNOrg)
    if orgGlycan:
        lines.append("| Rank | Organism | Sugar Sequence | Count |")
        lines.append("|:----:|----------|----------------|------:|")
        for rank, (org, seq, count) in enumerate(orgGlycan, 1):
            orgDisplay = org[:40] + "..." if len(org) > 40 else org
            lines.append(f"| {rank} | {orgDisplay} | `{seq}` | {count:,} |")
    else:
        lines.append("> 无物种来源数据可用\n")

    # Footer
    lines.append("\n---")
    lines.append("> Generated by `analyze_glycan_frequencies.py` | GlycoNP Pipeline\n")

    # Write
    reportContent = "\n".join(lines)
    with open(outputPath, "w", encoding="utf-8") as f:
        f.write(reportContent)

    print(f"  Report saved: {outputPath}")
    print(f"  Total lines: {len(lines)}")
    return outputPath


# =====================================================================
# Main
# =====================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="GlycoNP Glycomics Frequency Analysis")
    parser.add_argument("--input", type=str, default=None, help="Input pipeline CSV")
    parser.add_argument("--output", type=str, default=None, help="Output Markdown report")
    args = parser.parse_args()

    baseDir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

    if args.input:
        inputPath = args.input
    else:
        # 默认使用全量输出 (Default: full pipeline output)
        inputPath = os.path.join(baseDir, "reports", "GlycoNP_Pipeline_Full.csv")
        if not os.path.exists(inputPath):
            # 退回到 1000 条测试输出
            inputPath = os.path.join(baseDir, "reports", "GlycoNP_Pipeline_1000.csv")

    if args.output:
        outputPath = args.output
    else:
        outputPath = os.path.join(baseDir, "reports", "Glycomics_Summary.md")

    print("=" * 70)
    print("  GlycoNP Glycomics Frequency Analysis")
    print("=" * 70)
    print(f"  Input: {inputPath}")

    df = pd.read_csv(inputPath, low_memory=False, dtype=str, encoding="utf-8-sig")
    # 修正 bool 列 (Fix boolean columns from CSV string)
    for boolCol in ["Has_Nucleotide", "Has_Peptide"]:
        if boolCol in df.columns:
            df[boolCol] = df[boolCol].map({"True": True, "False": False, True: True, False: False}).fillna(False)

    print(f"  Loaded {len(df):,} rows")

    generateReport(df, outputPath)
    print(f"\n  Done!")
