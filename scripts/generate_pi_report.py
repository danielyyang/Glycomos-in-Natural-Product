"""
GlycoNP PI Report — 天然糖缀合物生物合成法则分析
(Biosynthetic Glycosylation Rules Report)

生成 4 个模块的图表 + JSON 统计数据:
  Module 1: Taxonomy × Sequence Heatmap
  Module 2: Aglycon → Sugar1 → Sugar2 Sankey
  Module 3: Modification Diversity Donut
  Module 4: Summary statistics for text report

使用方法 (Usage):
  python scripts/generate_pi_report.py
"""
import json
import os
import re
import sys
from collections import Counter, defaultdict

import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))


def cleanSuperclass(val: str) -> str:
    s = str(val).strip()
    if s in ("nan", "None", ""):
        return ""
    if "(Tanimoto=" in s:
        s = s[:s.index("(Tanimoto=")].strip()
    return s


def classifyKingdom(organisms: str) -> str:
    """根据物种名粗分界 (Rough kingdom classification)"""
    if not organisms or str(organisms) in ("nan", ""):
        return "Unknown"
    orgs = str(organisms).lower()
    # 细菌标记词
    bacteria = ["bacillus","streptomyces","pseudomonas","escherichia",
                "staphylococcus","mycobacterium","lactobacillus","clostridium",
                "salmonella","vibrio","enterococcus","actinomyces","nocardia",
                "corynebacterium","klebsiella","acinetobacter"]
    fungi = ["aspergillus","penicillium","fusarium","candida","saccharomyces",
             "trichoderma","neurospora","alternaria","mucor","rhizopus",
             "pleurotus","ganoderma","agaricus","trametes","boletus"]
    marine = ["sponge","coral","tunicate","ascidian","sea cucumber",
              "holothuria","strongylocentrotus","aplysina","halichondria",
              "axinella","ircinia","dysidea","haliclona"]
    if any(b in orgs for b in bacteria):
        return "Bacteria"
    if any(f in orgs for f in fungi):
        return "Fungi"
    if any(m in orgs for m in marine):
        return "Marine"
    return "Plantae"


def parseSequenceTokens(seq: str) -> list:
    """
    从 Sugar_Sequence 中提取有序的糖名列表。
    Extract ordered sugar name list from Sugar_Sequence.
    e.g., 'D-Glc-(?1-2)-D-GlcA' → ['D-Glc', 'D-GlcA']
    """
    if not seq or str(seq) in ("nan", "", "None"):
        return []
    # 清理预测后缀 (Clean prediction suffixes for display)
    cleanSeq = re.sub(r'_(?:predicted|iupac_parsed|biological_inferred|2d_inferred)', '', str(seq))
    tokens = re.findall(
        r'Neu5Ac|Neu5Gc|KDO|Non|Oct|Hept|'
        r'[DL]-[A-Z][a-z]+[A-Z]?[a-z]*(?:\([^)]*\))?|'
        r'Hex|dHex|Pen|HexA',
        cleanSeq)
    return tokens


def sequenceSignature(tokens: list, maxLen: int = 3) -> str:
    """
    将糖 token 列表转为摘要签名。
    Convert sugar token list to a short signature string.
    """
    if not tokens:
        return "N/A"
    truncated = tokens[:maxLen]
    return " → ".join(truncated)


def main():
    baseDir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    reportDir = os.path.join(baseDir, "reports")
    inputPath = os.path.join(reportDir, "GlycoNP_Pipeline_Full_Cleaned.csv")

    print("=" * 70)
    print("  GlycoNP PI Report — Biosynthetic Glycosylation Rules")
    print("  天然糖缀合物生物合成法则分析报告")
    print("=" * 70)

    df = pd.read_csv(inputPath, low_memory=False, encoding="utf-8-sig", dtype=str)
    # 过滤简单苷 (Filter simple glycosides)
    df = df[df["Is_Simple_Glycoside"] != "True"].copy()
    total = len(df)
    print(f"  Loaded: {total:,} compounds (excluded simple glycosides)")

    # 解析糖序列 (Parse sugar sequences)
    df["_tokens"] = df["Sugar_Sequence"].apply(parseSequenceTokens)
    df["_seq_sig"] = df["_tokens"].apply(sequenceSignature)
    df["_sugar1"] = df["_tokens"].apply(lambda t: t[0] if len(t) >= 1 else "N/A")
    df["_sugar2"] = df["_tokens"].apply(lambda t: t[1] if len(t) >= 2 else "—")
    df["_kingdom"] = df["organisms"].apply(classifyKingdom)
    df["_class"] = df["Superclass"].apply(cleanSuperclass)

    hasSugar = df["_sugar1"] != "N/A"
    sugarDf = df[hasSugar].copy()
    print(f"  Compounds with sugar(s): {len(sugarDf):,}")

    stats = {}

    # =================================================================
    # MODULE 1: Taxonomy × Sequence Heatmap
    # =================================================================
    print("\n  [Module 1] Taxonomy × Sequence Heatmap...")

    # Top 15 sequences
    seqCounts = sugarDf["_seq_sig"].value_counts()
    topSeqs = seqCounts.head(15).index.tolist()

    # Top 10 families (use LOTUS_Family first, fallback to kingdom)
    famCol = "LOTUS_Family"
    # For taxonomy dimension: use Kingdom for broader coverage
    kingdomCounts = sugarDf["_kingdom"].value_counts()
    topKingdoms = kingdomCounts.head(5).index.tolist()

    # Build heatmap: Kingdom × Sequence
    heatData = []
    for kingdom in topKingdoms:
        kDf = sugarDf[sugarDf["_kingdom"] == kingdom]
        kTotal = len(kDf)
        for seq in topSeqs:
            count = (kDf["_seq_sig"] == seq).sum()
            pct = round(count / kTotal * 100, 1) if kTotal > 0 else 0
            heatData.append({
                "Kingdom": kingdom,
                "Sequence": seq,
                "Pct": pct,
                "Count": count,
            })

    heatDf = pd.DataFrame(heatData)
    pivot = heatDf.pivot(index="Kingdom", columns="Sequence", values="Pct")
    pivot = pivot.reindex(columns=topSeqs)

    fig, ax = plt.subplots(figsize=(16, 6))
    sns.heatmap(pivot, annot=True, fmt=".1f", cmap="YlOrRd",
                linewidths=0.5, ax=ax, cbar_kws={"label": "% of Kingdom"})
    ax.set_title("Module 1: Taxonomic Glycosylation Patterns\n"
                 "Kingdom × Sugar Sequence Frequency (%)", fontsize=14, fontweight="bold")
    ax.set_xlabel("Sugar Sequence (Top 15)", fontsize=11)
    ax.set_ylabel("Kingdom", fontsize=11)
    plt.xticks(rotation=45, ha="right", fontsize=9)
    plt.tight_layout()
    heatmapPath = os.path.join(reportDir, "PI_Module1_Taxonomy_Heatmap.png")
    plt.savefig(heatmapPath, dpi=200)
    plt.close()
    print(f"    Saved: {heatmapPath}")

    # Kingdom-level stats for report
    stats["module1"] = {}
    for kingdom in topKingdoms:
        kDf = sugarDf[sugarDf["_kingdom"] == kingdom]
        top3 = kDf["_seq_sig"].value_counts().head(3)
        stats["module1"][kingdom] = {
            "total": len(kDf),
            "top3": [(seq, int(cnt), round(cnt/len(kDf)*100, 1))
                     for seq, cnt in top3.items()],
        }

    # =================================================================
    # MODULE 2: Aglycon → Sugar1 → Sugar2 Sankey
    # =================================================================
    print("\n  [Module 2] Aglycon → Sugar1 → Sugar2 Sankey...")

    targetClasses = ["Flavonoids", "Triterpenoids", "Steroids",
                     "Alkaloids", "Phenylpropanoids"]
    sankeyDf = sugarDf[sugarDf["_class"].isin(targetClasses)].copy()
    print(f"    Target classes: {len(sankeyDf):,} compounds")

    # Build Sankey data
    labels = []
    sources = []
    targets = []
    values = []
    colors = []

    # Aglycon → Sugar1 flows
    classColors = {
        "Flavonoids": "rgba(255,127,14,0.4)",
        "Triterpenoids": "rgba(44,160,44,0.4)",
        "Steroids": "rgba(148,103,189,0.4)",
        "Alkaloids": "rgba(214,39,40,0.4)",
        "Phenylpropanoids": "rgba(31,119,180,0.4)",
    }

    # Collect all unique labels
    allLabels = set()
    flowData = []
    for cls in targetClasses:
        clsDf = sankeyDf[sankeyDf["_class"] == cls]
        if len(clsDf) == 0:
            continue
        allLabels.add(cls)

        # Sugar1 distribution
        s1Counts = clsDf["_sugar1"].value_counts().head(5)
        for s1, cnt in s1Counts.items():
            s1Label = f"S1: {s1}"
            allLabels.add(s1Label)
            flowData.append((cls, s1Label, int(cnt), classColors.get(cls, "rgba(128,128,128,0.3)")))

            # Sugar2 distribution for this Sugar1
            s2Df = clsDf[(clsDf["_sugar1"] == s1) & (clsDf["_sugar2"] != "—")]
            if len(s2Df) > 0:
                s2Counts = s2Df["_sugar2"].value_counts().head(3)
                for s2, cnt2 in s2Counts.items():
                    s2Label = f"S2: {s2}"
                    allLabels.add(s2Label)
                    flowData.append((s1Label, s2Label, int(cnt2),
                                    classColors.get(cls, "rgba(128,128,128,0.3)")))

    labels = sorted(allLabels)
    labelIdx = {l: i for i, l in enumerate(labels)}

    for src, tgt, val, color in flowData:
        sources.append(labelIdx[src])
        targets.append(labelIdx[tgt])
        values.append(val)
        colors.append(color)

    fig = go.Figure(data=[go.Sankey(
        node=dict(
            pad=15,
            thickness=20,
            line=dict(color="black", width=0.5),
            label=labels,
            color=["rgba(100,100,100,0.8)"] * len(labels),
        ),
        link=dict(
            source=sources,
            target=targets,
            value=values,
            color=colors,
        )
    )])
    fig.update_layout(
        title_text="Module 2: Aglycon → Sugar1 → Sugar2 Flow<br>"
                   "<sub>Biosynthetic assembly paths for 5 major compound classes</sub>",
        font=dict(size=12, family="Segoe UI, Arial"),
        width=1400, height=800,
    )
    sankeyPath = os.path.join(reportDir, "PI_Module2_Aglycon_Sankey.html")
    fig.write_html(sankeyPath)
    print(f"    Saved: {sankeyPath}")

    # Stats for Module 2
    stats["module2"] = {}
    for cls in targetClasses:
        clsDf = sankeyDf[sankeyDf["_class"] == cls]
        if len(clsDf) == 0:
            continue
        s1Top = clsDf["_sugar1"].value_counts().head(3)
        s1TopSugar = s1Top.index[0] if len(s1Top) > 0 else "N/A"
        s1TopPct = round(s1Top.iloc[0] / len(clsDf) * 100, 1) if len(s1Top) > 0 else 0

        # 第二连接糖 (Sugar2)
        s2Df = clsDf[clsDf["_sugar2"] != "—"]
        s2Top = s2Df["_sugar2"].value_counts().head(3) if len(s2Df) > 0 else pd.Series(dtype=int)

        stats["module2"][cls] = {
            "total": len(clsDf),
            "s1_top": [(s, int(c), round(c/len(clsDf)*100,1)) for s,c in s1Top.items()],
            "s2_count": len(s2Df),
            "s2_pct": round(len(s2Df)/len(clsDf)*100, 1) if len(clsDf) > 0 else 0,
            "s2_top": [(s, int(c), round(c/len(s2Df)*100,1) if len(s2Df)>0 else 0)
                       for s,c in s2Top.head(3).items()],
        }

    # =================================================================
    # MODULE 3: Modification Diversity Donut
    # =================================================================
    print("\n  [Module 3] Modification Diversity...")

    # 判断是否有修饰 (Has modification?)
    # Glycan_Modifications 格式: D-Glc_1(*S) → *S = 有标记
    def hasRealModification(mods: str) -> bool:
        if not mods or str(mods) in ("nan", ""):
            return False
        # 有 Acetylated, Methylated, Sulfated, Phosphorylated 等
        modStr = str(mods).lower()
        realMods = ["acetyl", "methyl", "sulfat", "phosph", "carbamoyl",
                    "hydroxyl", "glycosyl", "caffeoyl", "coumaroyl",
                    "feruloyl", "malonyl", "sinapoyl"]
        return any(m in modStr for m in realMods)

    sugarDf["_has_mod"] = sugarDf["Glycan_Modifications"].apply(hasRealModification)
    modRate = sugarDf["_has_mod"].sum() / len(sugarDf) * 100

    # 按大类统计修饰率
    modByClass = []
    topClasses = sugarDf["_class"].value_counts().head(10).index.tolist()
    topClasses = [c for c in topClasses if c and c != "Unclassified"]

    for cls in topClasses:
        clsDf = sugarDf[sugarDf["_class"] == cls]
        if len(clsDf) < 20:
            continue
        rate = clsDf["_has_mod"].sum() / len(clsDf) * 100
        modByClass.append({"Class": cls, "Modified (%)": round(rate, 1),
                           "Bare (%)": round(100-rate, 1), "Total": len(clsDf)})

    modClassDf = pd.DataFrame(modByClass)

    if len(modClassDf) > 0:
        # Stacked bar
        fig, ax = plt.subplots(figsize=(12, 6))
        x = range(len(modClassDf))
        ax.barh(modClassDf["Class"], modClassDf["Modified (%)"],
                color="#e74c3c", label="Modified", alpha=0.85)
        ax.barh(modClassDf["Class"], modClassDf["Bare (%)"],
                left=modClassDf["Modified (%)"],
                color="#3498db", label="Bare (unmodified)", alpha=0.85)
        ax.set_xlabel("Percentage (%)", fontsize=11)
        ax.set_title("Module 3: Sugar Modification Rate by Compound Class\n"
                    "(修饰基团覆盖率: Acetyl, Methyl, Sulfate, Phosphate, etc.)",
                    fontsize=13, fontweight="bold")
        ax.legend(loc="lower right")
        for i, row in modClassDf.iterrows():
            ax.text(row["Modified (%)"] / 2, i, f'{row["Modified (%)"]}%',
                    ha="center", va="center", fontsize=9, color="white", fontweight="bold")
        plt.tight_layout()
        donutPath = os.path.join(reportDir, "PI_Module3_Modification_Bar.png")
        plt.savefig(donutPath, dpi=200)
        plt.close()
        print(f"    Saved: {donutPath}")

    # Overall donut
    fig, ax = plt.subplots(figsize=(7, 7))
    modCount = int(sugarDf["_has_mod"].sum())
    bareCount = len(sugarDf) - modCount
    ax.pie([modCount, bareCount],
           labels=[f"Modified\n{modCount:,} ({modRate:.1f}%)",
                   f"Bare\n{bareCount:,} ({100-modRate:.1f}%)"],
           colors=["#e74c3c", "#3498db"],
           autopct="", startangle=140,
           wedgeprops=dict(width=0.4, edgecolor="white"))
    ax.set_title("Overall Sugar Modification Rate\n(整体糖链修饰覆盖率)",
                fontsize=13, fontweight="bold")
    donutOverall = os.path.join(reportDir, "PI_Module3_Donut_Overall.png")
    plt.savefig(donutOverall, dpi=200)
    plt.close()
    print(f"    Saved: {donutOverall}")

    stats["module3"] = {
        "overall_mod_rate": round(modRate, 1),
        "modified_count": modCount,
        "bare_count": bareCount,
        "by_class": modByClass,
    }

    # =================================================================
    # MODULE 4: Summary Statistics
    # =================================================================
    print("\n  [Module 4] Summary Statistics...")

    # 高频特异性组合
    comboDf = sugarDf[sugarDf["_class"].isin(targetClasses)].copy()
    comboDf["_combo"] = comboDf["_class"] + " + " + comboDf["_seq_sig"]
    comboTop = comboDf["_combo"].value_counts().head(10)
    stats["module4"] = {
        "top_combos": [(combo, int(cnt)) for combo, cnt in comboTop.items()],
        "total_with_sugar": len(sugarDf),
        "multi_sugar_pct": round((sugarDf["_sugar2"] != "—").sum() / len(sugarDf) * 100, 1),
    }

    # 保存 JSON
    jsonPath = os.path.join(reportDir, "PI_Report_Stats.json")
    with open(jsonPath, "w", encoding="utf-8") as f:
        json.dump(stats, f, ensure_ascii=False, indent=2)
    print(f"    Saved: {jsonPath}")

    print(f"\n{'='*70}")
    print(f"  PI Report Generation Complete!")
    print(f"  Charts: 4 files")
    print(f"  Stats JSON: {jsonPath}")
    print(f"{'='*70}")

    # Print summary for report compilation
    print("\n  === REPORT DATA SUMMARY ===")
    print(f"\n  Module 1 - Kingdom top sequences:")
    for k, v in stats["module1"].items():
        print(f"    {k} ({v['total']:,}):")
        for seq, cnt, pct in v["top3"]:
            print(f"      {seq:30s} {cnt:>5,} ({pct}%)")

    print(f"\n  Module 2 - Aglycon sugar preferences:")
    for cls, v in stats["module2"].items():
        print(f"    {cls} ({v['total']:,}):")
        print(f"      Sugar1 top: {v['s1_top']}")
        print(f"      {v['s2_pct']}% have Sugar2, top: {v['s2_top'][:2]}")

    print(f"\n  Module 3 - Modification rate: {stats['module3']['overall_mod_rate']}%")

    print(f"\n  Module 4 - Top combos:")
    for combo, cnt in stats["module4"]["top_combos"][:5]:
        print(f"    {combo:50s} {cnt:>5,}")


if __name__ == "__main__":
    main()
