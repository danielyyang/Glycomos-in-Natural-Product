"""
GlycoNP 高级可视化 — 主任汇报级别图表
GlycoNP Executive Plots — PI-level Presentation Charts

三大图表 (Three Charts):
  1. 桑基图 (Sankey): Superclass → Modification → Sugar_Sequence — 糖链去哪儿了？
  2. 热力图 (Heatmap): Sugar_Sequence × Superclass — 大自然的黄金搭档
  3. 分组柱状图 (Grouped Bar): Kingdom × Modification — 修饰的跨界印记

使用方法 (Usage):
  python scripts/generate_executive_plots.py [--input PATH]
"""
import argparse
import os
import sys
from collections import Counter

import numpy as np
import pandas as pd

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))


# =====================================================================
# 工具函数 (Utility Functions)
# =====================================================================

def cleanSuperclass(val: str) -> str:
    """
    清理 Superclass 名称: 去掉 Tanimoto 后缀。
    Clean Superclass name: remove Tanimoto suffix.
    "Triterpenoids (Tanimoto=1.000)" → "Triterpenoids"
    """
    if not val or val == "nan":
        return "Unclassified"
    s = str(val).strip()
    if "(Tanimoto=" in s:
        s = s[:s.index("(Tanimoto=")].strip()
    return s if s else "Unclassified"


def inferKingdom(organisms: str) -> str:
    """
    从物种名推断生物界 (Kingdom)。
    Infer Kingdom from organism names.

    设计原则: 先精确匹配关键词, 再用常见属名推断。
    """
    if not organisms or str(organisms) == "nan":
        return "Unknown"

    text = str(organisms).lower()

    # 细菌标志 (Bacteria markers)
    bacteriaGenera = [
        "streptomyces", "bacillus", "pseudomonas", "escherichia",
        "staphylococcus", "salmonella", "mycobacterium", "clostridium",
        "lactobacillus", "enterococcus", "vibrio", "actinomyces",
        "nocardia", "micromonospora", "saccharopolyspora",
    ]
    if any(g in text for g in bacteriaGenera):
        return "Bacteria"

    # 真菌标志 (Fungi markers)
    fungiGenera = [
        "aspergillus", "penicillium", "fusarium", "candida",
        "saccharomyces", "trichoderma", "agaricus", "ganoderma",
        "amanita", "botrytis", "cladosporium",
    ]
    if any(g in text for g in fungiGenera):
        return "Fungi"

    # 海洋生物标志 (Marine markers)
    marineIndicators = [
        "sponge", "coral", "sea", "marine", "halichondria",
        "aplysina", "haliclona", "ircinia",
    ]
    if any(m in text for m in marineIndicators):
        return "Marine"

    # 默认植物 (Default: Plantae — most COCONUT entries are plants)
    return "Plantae"


def parseModifications(modStr: str) -> list:
    """
    解析 Glycan_Modifications 字符串, 返回修饰名列表。
    Parse modification string → list of modification names.
    "O-Ac: 2, NAc: 1" → ["O-Ac", "O-Ac", "NAc"]
    """
    if not modStr or str(modStr) in ("nan", ""):
        return []

    mods = []
    for part in str(modStr).split(","):
        part = part.strip()
        if ":" in part:
            name = part.split(":")[0].strip()
            try:
                count = int(part.split(":")[1].strip())
            except ValueError:
                count = 1
            mods.extend([name] * count)
    return mods


# =====================================================================
# 图表 1: 桑基图 (Sankey Diagram)
# =====================================================================

def plotSankey(df: pd.DataFrame, outputPath: str):
    """
    绘制 Superclass → Modification → Sugar_Sequence 桑基图。
    Sankey: Superclass → Modification → Sugar_Sequence.

    "糖链去哪儿了？" — Where Do the Sugars Go?
    """
    import plotly.graph_objects as go

    print("\n[Chart 1] Sankey Diagram: Where Do the Sugars Go?")

    # 准备数据: 取 Top 10 Superclass
    df["_class"] = df["Superclass"].apply(cleanSuperclass)
    topClasses = df["_class"].value_counts().head(10).index.tolist()
    subDf = df[df["_class"].isin(topClasses)].copy()

    # 准备修饰标签 (Extract primary modification)
    def primaryMod(modStr):
        mods = parseModifications(modStr)
        if not mods:
            return "Unmodified"
        # 取出现最多的修饰
        return Counter(mods).most_common(1)[0][0]

    subDf["_mod"] = subDf["Glycan_Modifications"].apply(primaryMod)

    # 准备糖序列 (Clean sugar sequence)
    subDf["_seq"] = subDf["Sugar_Sequence"].fillna("Unknown").astype(str)
    subDf.loc[subDf["_seq"].isin(["", "nan"]), "_seq"] = "Unknown"
    # 取全局 Top 8 序列, 其余归 "Other"
    topSeqs = subDf["_seq"].value_counts().head(8).index.tolist()
    subDf.loc[~subDf["_seq"].isin(topSeqs), "_seq"] = "Other"

    # 构建节点 (Build nodes)
    classNodes = sorted(topClasses)
    modTypes = sorted(subDf["_mod"].unique())
    seqNodes = sorted(subDf["_seq"].unique())

    allNodes = classNodes + modTypes + seqNodes
    nodeIdx = {name: i for i, name in enumerate(allNodes)}

    # 颜色方案 (Color scheme)
    # 大类: 蓝色调, 修饰: 黄色调, 糖序列: 红色调
    nClass = len(classNodes)
    nMod = len(modTypes)
    nSeq = len(seqNodes)

    classColors = [f"hsla({int(210 + i*15)}, 70%, 55%, 0.8)" for i in range(nClass)]
    modColors = [f"hsla({int(40 + i*12)}, 85%, 55%, 0.8)" for i in range(nMod)]
    seqColors = [f"hsla({int(0 + i*15)}, 70%, 60%, 0.8)" for i in range(nSeq)]
    nodeColors = classColors + modColors + seqColors

    # 构建连接 (Build links): Class → Mod
    linkSrc, linkTgt, linkVal, linkColor = [], [], [], []

    # Layer 1: Superclass → Modification
    for cls in classNodes:
        clsDf = subDf[subDf["_class"] == cls]
        modCounts = clsDf["_mod"].value_counts()
        for mod, count in modCounts.items():
            if count > 0:
                linkSrc.append(nodeIdx[cls])
                linkTgt.append(nodeIdx[mod])
                linkVal.append(count)
                linkColor.append(f"hsla(210, 50%, 70%, 0.3)")

    # Layer 2: Modification → Sugar_Sequence
    for mod in modTypes:
        modDf = subDf[subDf["_mod"] == mod]
        seqCounts = modDf["_seq"].value_counts()
        for seq, count in seqCounts.items():
            if count > 0:
                linkSrc.append(nodeIdx[mod])
                linkTgt.append(nodeIdx[seq])
                linkVal.append(count)
                linkColor.append(f"hsla(40, 60%, 70%, 0.3)")

    fig = go.Figure(go.Sankey(
        arrangement="snap",
        node=dict(
            pad=18,
            thickness=25,
            line=dict(color="white", width=1),
            label=allNodes,
            color=nodeColors,
        ),
        link=dict(
            source=linkSrc,
            target=linkTgt,
            value=linkVal,
            color=linkColor,
        ),
    ))

    fig.update_layout(
        title=dict(
            text="糖链去哪儿了？ — Where Do the Sugars Go?<br>"
                 "<sub>Superclass → Primary Modification → Sugar Sequence</sub>",
            font=dict(size=18),
        ),
        font=dict(size=12, family="Segoe UI, Arial"),
        width=1200, height=700,
        paper_bgcolor="#fafafa",
    )

    fig.write_html(outputPath)
    print(f"  Saved: {outputPath}")
    return outputPath


# =====================================================================
# 图表 2: 热力图 (Heatmap)
# =====================================================================

def plotHeatmap(df: pd.DataFrame, outputPath: str):
    """
    绘制 Sugar_Sequence × Superclass 交叉热力图。
    Heatmap: Sugar_Sequence × Superclass (log-scaled).

    "大自然的黄金搭档" — Nature's Golden Pairs.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns

    print("\n[Chart 2] Heatmap: Nature's Golden Pairs")

    df["_class"] = df["Superclass"].apply(cleanSuperclass)
    df["_seq"] = df["Sugar_Sequence"].fillna("Unknown").astype(str)
    df.loc[df["_seq"].isin(["", "nan"]), "_seq"] = "Unknown"

    # Top 20 sequences, Top 15 classes
    topSeqs = df["_seq"].value_counts().head(20).index.tolist()
    topClasses = df["_class"].value_counts().head(15).index.tolist()

    subDf = df[df["_seq"].isin(topSeqs) & df["_class"].isin(topClasses)]

    # 构建交叉表 (Pivot table)
    pivot = pd.crosstab(subDf["_class"], subDf["_seq"])
    # 按行总数排序 (Sort by row totals)
    pivot = pivot.loc[pivot.sum(axis=1).sort_values(ascending=False).index]
    # 按列总数排序 (Sort by column totals)
    pivot = pivot[pivot.sum(axis=0).sort_values(ascending=False).index]

    # Log 缩放 (Log scale: log2(x+1))
    pivotLog = np.log2(pivot + 1)

    # 绘图 (Draw)
    fig, ax = plt.subplots(figsize=(16, 10))

    sns.heatmap(
        pivotLog,
        annot=pivot.values,  # 原始频次作为标注
        fmt="d",
        cmap="YlOrRd",
        linewidths=0.5,
        linecolor="#f0f0f0",
        cbar_kws={"label": "log₂(count + 1)", "shrink": 0.7},
        ax=ax,
        annot_kws={"size": 8},
    )

    ax.set_title(
        "Nature's Golden Pairs\n"
        "Sugar Sequence × Superclass Frequency (log2-scaled)",
        fontsize=16, fontweight="bold", pad=20,
    )
    ax.set_xlabel("Sugar Sequence", fontsize=12, labelpad=10)
    ax.set_ylabel("Superclass", fontsize=12, labelpad=10)

    plt.xticks(rotation=45, ha="right", fontsize=9)
    plt.yticks(fontsize=10)
    plt.tight_layout()

    fig.savefig(outputPath, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved: {outputPath}")
    return outputPath


# =====================================================================
# 图表 3: 分组柱状图 (Grouped Bar Chart)
# =====================================================================

def plotKingdomMods(df: pd.DataFrame, outputPath: str):
    """
    绘制 Kingdom × Modification 分组柱状图。
    Grouped Bar: Kingdom × Modification.

    "修饰的跨界印记" — Cross-Kingdom Modification Fingerprint.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    print("\n[Chart 3] Grouped Bar: Cross-Kingdom Modification Fingerprint")

    # 推断 Kingdom (Infer Kingdom)
    df["_kingdom"] = df["organisms"].apply(inferKingdom)
    kingdomCounts = df["_kingdom"].value_counts()
    print(f"  Kingdom distribution: {dict(kingdomCounts)}")

    # 只保留有足够数据的 Kingdom
    validKingdoms = [k for k in ["Plantae", "Bacteria", "Fungi", "Marine"]
                     if k in kingdomCounts.index and kingdomCounts[k] >= 5]

    if not validKingdoms:
        validKingdoms = kingdomCounts.head(4).index.tolist()

    targetMods = ["O-Me", "O-Ac", "NAc", "Sulfate", "Phosphate", "COOH", "NH2"]

    # 统计每个 Kingdom 中含各修饰的化合物比例
    dataRows = []
    for kingdom in validKingdoms:
        kDf = df[df["_kingdom"] == kingdom]
        kTotal = len(kDf)
        for mod in targetMods:
            # 含有该修饰的化合物数
            hasMod = sum(1 for m in kDf["Glycan_Modifications"].fillna("")
                         if mod in str(m))
            pct = hasMod / kTotal * 100 if kTotal > 0 else 0
            dataRows.append({"Kingdom": kingdom, "Modification": mod,
                             "Percentage": pct, "Count": hasMod})

    plotDf = pd.DataFrame(dataRows)

    # 绘图 (Draw)
    fig, ax = plt.subplots(figsize=(14, 7))

    kingdoms = validKingdoms
    nKingdoms = len(kingdoms)
    nMods = len(targetMods)
    barWidth = 0.8 / nKingdoms
    x = np.arange(nMods)

    # 颜色方案 — 每个界一种颜色
    kingdomColors = {
        "Plantae": "#2ecc71",
        "Bacteria": "#e74c3c",
        "Fungi": "#9b59b6",
        "Marine": "#3498db",
        "Unknown": "#95a5a6",
    }

    for i, kingdom in enumerate(kingdoms):
        kData = plotDf[plotDf["Kingdom"] == kingdom]
        vals = [kData[kData["Modification"] == m]["Percentage"].values[0]
                if len(kData[kData["Modification"] == m]) > 0 else 0
                for m in targetMods]
        color = kingdomColors.get(kingdom, "#7f8c8d")
        bars = ax.bar(x + i * barWidth, vals, barWidth,
                      label=kingdom, color=color, edgecolor="white",
                      linewidth=0.5, alpha=0.85)

        # 数值标签 (Value labels)
        for bar, val in zip(bars, vals):
            if val > 0.5:
                ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.3,
                        f"{val:.1f}%", ha="center", va="bottom",
                        fontsize=7, fontweight="bold", color="#333")

    ax.set_xlabel("Modification Type", fontsize=13, labelpad=10)
    ax.set_ylabel("% of Compounds with Modification", fontsize=13, labelpad=10)
    ax.set_title(
        "Cross-Kingdom Modification Fingerprint\n"
        "Proportion of compounds carrying each modification, by biological kingdom",
        fontsize=15, fontweight="bold", pad=15,
    )
    ax.set_xticks(x + barWidth * (nKingdoms - 1) / 2)
    ax.set_xticklabels(targetMods, fontsize=11)
    ax.legend(title="Kingdom", fontsize=11, title_fontsize=12,
              loc="upper right", framealpha=0.9)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", alpha=0.3, linestyle="--")

    plt.tight_layout()
    fig.savefig(outputPath, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved: {outputPath}")
    return outputPath


# =====================================================================
# Main
# =====================================================================

def main():
    parser = argparse.ArgumentParser(description="GlycoNP Executive Plots")
    parser.add_argument("--input", type=str, default=None, help="Input CSV")
    args = parser.parse_args()

    baseDir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    reportDir = os.path.join(baseDir, "reports")

    if args.input:
        inputPath = args.input
    else:
        # 优先使用清洗后数据 (Prefer cleaned data)
        inputPath = os.path.join(reportDir, "GlycoNP_Pipeline_Full_Cleaned.csv")
        if not os.path.exists(inputPath):
            inputPath = os.path.join(reportDir, "GlycoNP_Pipeline_Full.csv")

    print("=" * 70)
    print("  GlycoNP Executive Plots")
    print("=" * 70)
    print(f"  Input: {inputPath}")

    df = pd.read_csv(inputPath, low_memory=False, dtype=str, encoding="utf-8-sig")
    # 修正 bool 列
    for boolCol in ["Has_Nucleotide", "Has_Peptide", "Is_Simple_Glycoside"]:
        if boolCol in df.columns:
            df[boolCol] = df[boolCol].map(
                {"True": True, "False": False, True: True, False: False}
            ).fillna(False)

    # 过滤假苷元 (Filter out simple glycosides)
    nBefore = len(df)
    if "Is_Simple_Glycoside" in df.columns:
        df = df[df["Is_Simple_Glycoside"] != True].copy()
    print(f"  Loaded {nBefore:,} → filtered to {len(df):,} rows (excluded {nBefore-len(df):,} simple glycosides)")

    # Chart 1: Sankey
    sankeyPath = os.path.join(reportDir, "Chart_Sankey_SugarFlow.html")
    plotSankey(df, sankeyPath)

    # Chart 2: Heatmap
    heatmapPath = os.path.join(reportDir, "Chart_Heatmap_GoldenPairs.png")
    plotHeatmap(df, heatmapPath)

    # Chart 3: Grouped Bar
    barPath = os.path.join(reportDir, "Chart_Bar_CrossKingdom.png")
    plotKingdomMods(df, barPath)

    print(f"\n{'='*70}")
    print(f"  All 3 charts generated!")
    print(f"  1. {sankeyPath}")
    print(f"  2. {heatmapPath}")
    print(f"  3. {barPath}")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
