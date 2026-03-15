"""
GlycoNP 全大类三维配方图谱 — 骨架×糖链×物种 可视化
(Full Superclass Scaffold × Sugar × Taxonomy Visualization)

自动发现所有化合物大类, 逐一生成:
  1. Sunburst: Murcko_Scaffold → Sugar_Sequence → LOTUS_Family
  2. Bubble Chart: Sugar × Family grid, colored by scaffold

使用方法 (Usage):
  python scripts/plot_scaffold_taxonomy_map.py
"""
import json
import os
import re
import sys
from collections import Counter, OrderedDict

import pandas as pd
import plotly.express as px

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

# 最低化合物数阈值 (Minimum compounds threshold for generating charts)
MIN_COMPOUNDS = 30
TOP_SCAFFOLDS = 10


def cleanSuperclass(val: str) -> str:
    """清理 Superclass 名称 (strip Tanimoto tags)."""
    s = str(val).strip()
    if s in ("nan", "None", ""):
        return ""
    if "(Tanimoto=" in s:
        s = s[:s.index("(Tanimoto=")].strip()
    return s


def shortenScaffold(smiles: str, maxLen: int = 35) -> str:
    """缩短骨架 SMILES 显示"""
    s = str(smiles)
    return s[:maxLen] + "…" if len(s) > maxLen else s


def extractFirstSugar(seq: str) -> str:
    """从糖序列中提取最主要的糖名"""
    if not seq or str(seq) in ("nan", "", "None"):
        return "N/A"
    tokens = re.findall(
        r'Neu5Ac|Neu5Gc|KDO|Non|Oct|Hept|'
        r'[DL]-[A-Z][a-z]+[A-Z]?[a-z]*(?:\([^)]*\))?|'
        r'Hex|dHex|Pen|HexA',
        str(seq))
    if tokens:
        return tokens[0]
    return str(seq)[:20]


def processSingleClass(
    classDf: pd.DataFrame,
    className: str,
    reportDir: str,
) -> dict:
    """
    处理单个大类: 生成 Sunburst + Bubble, 提取 Top-1 骨架糖偏好。
    Process single class: generate charts + extract scaffold insights.
    Returns dict with analysis results.
    """
    result = {
        "class": className,
        "totalCompounds": len(classDf),
        "chartGenerated": False,
        "top1Scaffold": None,
        "top1ScaffoldSmiles": None,
        "top1Count": 0,
        "top1Pct": 0,
        "sugarPrefs": [],
    }

    if len(classDf) < MIN_COMPOUNDS:
        return result

    # 有效性过滤: 需要骨架和糖链 (Need scaffold + sugar)
    validDf = classDf[
        classDf["Murcko_Scaffold"].notna()
        & (classDf["Murcko_Scaffold"] != "")
        & (classDf["Murcko_Scaffold"] != "nan")
        & classDf["Sugar_Sequence"].notna()
        & (classDf["Sugar_Sequence"] != "")
        & (classDf["Sugar_Sequence"] != "nan")
    ].copy()

    if len(validDf) < 10:
        return result

    # Top N 骨架
    scaffoldCounts = validDf["Murcko_Scaffold"].value_counts()
    topScaffolds = scaffoldCounts.head(TOP_SCAFFOLDS).index.tolist()
    chartDf = validDf[validDf["Murcko_Scaffold"].isin(topScaffolds)].copy()

    chartDf["Scaffold_Short"] = chartDf["Murcko_Scaffold"].apply(
        lambda x: shortenScaffold(x, 35))
    chartDf["Primary_Sugar"] = chartDf["Sugar_Sequence"].apply(extractFirstSugar)
    chartDf["Family"] = chartDf["LOTUS_Family"].fillna("Unknown").apply(
        lambda x: str(x)[:30] if str(x) not in ("nan", "None", "") else "Unknown"
    )

    # ===== Sunburst =====
    aggSun = (
        chartDf.groupby(["Scaffold_Short", "Primary_Sugar", "Family"])
        .size()
        .reset_index(name="Count")
    )
    aggSun = aggSun[aggSun["Count"] >= 2]

    # 文件名安全处理 (Sanitize class name for filename)
    safeClassName = re.sub(r'[^\w\-]', '_', className)

    if len(aggSun) > 0:
        fig = px.sunburst(
            aggSun,
            path=["Scaffold_Short", "Primary_Sugar", "Family"],
            values="Count",
            title=f"{className} — Scaffold × Sugar × Family Sunburst<br>"
                  f"<sub>Top {TOP_SCAFFOLDS} Murcko scaffolds | "
                  f"{len(chartDf):,} compounds</sub>",
            color="Count",
            color_continuous_scale="YlOrRd",
        )
        fig.update_layout(
            width=1200, height=900,
            font=dict(family="Segoe UI, Arial", size=13),
            margin=dict(t=80, l=10, r=10, b=10),
        )
        sunPath = os.path.join(reportDir,
                               f"Sunburst_Scaffold_{safeClassName}.html")
        fig.write_html(sunPath)
        result["chartGenerated"] = True

    # ===== Bubble =====
    topSugars = chartDf["Primary_Sugar"].value_counts().head(15).index.tolist()
    topFamilies = chartDf["Family"].value_counts().head(15).index.tolist()
    bubbleDf = chartDf[
        chartDf["Primary_Sugar"].isin(topSugars)
        & chartDf["Family"].isin(topFamilies)
    ].copy()

    aggBub = (
        bubbleDf.groupby(["Primary_Sugar", "Family", "Scaffold_Short"])
        .size()
        .reset_index(name="Count")
    )
    aggBub = aggBub[aggBub["Count"] >= 2]

    if len(aggBub) > 0:
        fig2 = px.scatter(
            aggBub,
            x="Primary_Sugar",
            y="Family",
            size="Count",
            color="Scaffold_Short",
            hover_data={"Count": True, "Scaffold_Short": True},
            size_max=40,
            title=f"{className} — Sugar × Family × Scaffold Bubble Map<br>"
                  f"<sub>Bubble size = compound count | "
                  f"Color = Murcko scaffold</sub>",
            labels={
                "Primary_Sugar": "Sugar Type",
                "Family": "Taxonomic Family",
                "Scaffold_Short": "Scaffold",
            },
        )
        fig2.update_layout(
            width=1400, height=900,
            font=dict(family="Segoe UI, Arial", size=12),
            xaxis=dict(tickangle=-45),
            legend=dict(font=dict(size=9), title="Scaffold"),
            margin=dict(b=120),
        )
        bubPath = os.path.join(reportDir,
                                f"BubbleMap_{safeClassName}.html")
        fig2.write_html(bubPath)
        result["chartGenerated"] = True

    # ===== Top-1 骨架分析 =====
    top1Smiles = topScaffolds[0]
    top1Df = chartDf[chartDf["Murcko_Scaffold"] == top1Smiles]
    sugarCounts = top1Df["Primary_Sugar"].value_counts()

    result["top1ScaffoldSmiles"] = top1Smiles
    result["top1Scaffold"] = shortenScaffold(top1Smiles, 60)
    result["top1Count"] = len(top1Df)
    result["top1Pct"] = round(len(top1Df) / len(chartDf) * 100, 1) \
        if len(chartDf) > 0 else 0

    for sugar, cnt in sugarCounts.head(3).items():
        pct = round(cnt / len(top1Df) * 100, 1) if len(top1Df) > 0 else 0
        result["sugarPrefs"].append({
            "sugar": sugar,
            "count": int(cnt),
            "pct": pct,
        })

    return result


def main():
    baseDir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    reportDir = os.path.join(baseDir, "reports")
    inputPath = os.path.join(reportDir, "GlycoNP_Pipeline_Full_Cleaned.csv")

    print("=" * 70)
    print("  GlycoNP Full Superclass Scaffold × Sugar × Taxonomy Map")
    print("  全大类骨架-糖链-物种三维配方图谱")
    print("=" * 70)

    df = pd.read_csv(inputPath, low_memory=False, dtype=str,
                     encoding="utf-8-sig")
    total = len(df)
    print(f"  Loaded: {total:,} rows")

    # 清理分类
    df["_class"] = df["Superclass"].apply(cleanSuperclass)

    # 动态发现所有大类 (Dynamic class discovery)
    allClasses = (
        df["_class"]
        .loc[df["_class"] != ""]
        .value_counts()
    )
    print(f"  Unique classes: {len(allClasses)}")

    # 排除 Unclassified 和极小类 (Exclude unclassified + tiny)
    targetClasses = [
        cls for cls, cnt in allClasses.items()
        if cnt >= MIN_COMPOUNDS and cls != "Unclassified"
    ]
    print(f"  Target classes (≥{MIN_COMPOUNDS} compounds, excl. Unclassified): "
          f"{len(targetClasses)}")

    results = []
    filesGenerated = 0

    for i, className in enumerate(targetClasses, 1):
        cnt = allClasses[className]
        print(f"\n  [{i}/{len(targetClasses)}] {className} ({cnt:,} compounds)")
        classDf = df[df["_class"] == className]

        result = processSingleClass(classDf, className, reportDir)
        results.append(result)

        if result["chartGenerated"]:
            filesGenerated += 2  # Sunburst + Bubble

    # ===== 输出 JSON 汇总 (Save JSON summary for report) =====
    jsonPath = os.path.join(reportDir, "scaffold_taxonomy_results.json")
    with open(jsonPath, "w", encoding="utf-8") as f:
        json.dump(results, f, ensure_ascii=False, indent=2)

    print(f"\n{'='*70}")
    print(f"  Complete!")
    print(f"  Classes processed: {len(targetClasses)}")
    print(f"  HTML files generated: {filesGenerated}")
    print(f"  Results JSON: {jsonPath}")
    print(f"{'='*70}")

    # ===== 终端汇报 (Terminal Report) =====
    print(f"\n{'='*70}")
    print(f"  TOP-1 SCAFFOLD SUGAR PREFERENCES — ALL CLASSES")
    print(f"{'='*70}")
    for r in results:
        if not r["chartGenerated"]:
            continue
        print(f"\n  {r['class']}:")
        print(f"    Top-1 scaffold: {r['top1Scaffold']}")
        print(f"    Compounds: {r['top1Count']} ({r['top1Pct']}%)")
        for j, sp in enumerate(r["sugarPrefs"], 1):
            medal = ["🥇", "🥈", "🥉"][j-1]
            print(f"    {medal} {sp['sugar']:20s} {sp['count']:>4} ({sp['pct']}%)")


if __name__ == "__main__":
    main()
