#!/usr/bin/env python3
"""
==========================================================================
  [EN] GlycoNP Saponin Database — Bioactivity-Oriented Analytics (F-Series)
       Charts cross-referencing glycan structure with pharmacological data:
       chain length vs bioactivity, modification vs pChEMBL, sugar vs target,
       Lipinski RO5 compliance, and QED drug-likeness.

  [CN] GlycoNP 皂苷子库 — 药学/生物活性导向分析图谱 (F 系列)
       糖链结构与药理学数据的交叉分析: 链长 × 生物活性, 修饰 × pChEMBL,
       糖类偏好 × 靶点蛋白, Lipinski RO5 合规性, QED 类药性。

  Input:  reports/GlycoNP_Saponin_DB.csv
  Output: reports/saponin_figures/F*.html + F*.png
  [TEST DATA ONLY]
==========================================================================
"""
import os
import re
import sys
import warnings
from pathlib import Path
from collections import Counter
from typing import List, Dict

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

warnings.filterwarnings("ignore")

# ── 全局变量: 由 _loadModule 注入或独立运行时自行定义 ─────────────────
try:
    _ = LAYOUT_DEFAULTS  # type: ignore
except NameError:
    BASE_DIR = Path(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
    CSV_PATH = BASE_DIR / "reports" / "GlycoNP_Saponin_DB.csv"
    OUT_DIR = BASE_DIR / "reports" / "saponin_figures"
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    LAYOUT_DEFAULTS = dict(
        template="plotly_white",
        font=dict(family="Arial", size=13, color="#222"),
        paper_bgcolor="white", plot_bgcolor="white",
        margin=dict(l=80, r=40, t=80, b=80),
        title_font=dict(size=17, family="Arial", color="#222"),
    )
    HEATMAP_COLORSCALE = "Viridis"
    QUAL_COLORS = px.colors.qualitative.Set2 + px.colors.qualitative.Pastel1

    def saveHtml(fig, filename: str, height: int = 650, width: int = 1100):
        """Save Plotly figure as HTML + PNG."""
        curMargin = fig.layout.margin
        customL, customR = getattr(curMargin, "l", None), getattr(curMargin, "r", None)
        customT, customB = getattr(curMargin, "t", None), getattr(curMargin, "b", None)
        fig.update_layout(**LAYOUT_DEFAULTS)
        fig.update_layout(margin=dict(
            l=customL if customL is not None else 80,
            r=customR if customR is not None else 40,
            t=customT if customT is not None else 80,
            b=customB if customB is not None else 80,
        ))
        if height: fig.update_layout(height=height)
        if width: fig.update_layout(width=width)
        outPath = OUT_DIR / filename
        fig.write_html(str(outPath), include_plotlyjs="cdn")
        print(f"  -> Saved: {outPath.name}")
        pngPath = OUT_DIR / filename.replace(".html", ".png")
        try:
            curMargin = fig.layout.margin
            fig.update_layout(margin=dict(
                l=max(getattr(curMargin, "l", 80) or 80, 100),
                r=max(getattr(curMargin, "r", 40) or 40, 60),
                t=max(getattr(curMargin, "t", 80) or 80, 90),
                b=max(getattr(curMargin, "b", 80) or 80, 100),
            ))
            fig.write_image(str(pngPath), format="png", scale=3,
                            width=width, height=height, engine="kaleido")
            print(f"  -> Saved: {pngPath.name}")
        except Exception as e:
            print(f"  [WARN] PNG export failed ({e}).")

    def extractMonosaccharideList(seq: str) -> List[str]:
        if pd.isna(seq) or not seq: return []
        s = re.sub(r"\([ab]\d-\d\)", " ", str(seq))
        s = s.replace(";", " ").replace("[", " ").replace("]", " ")
        return re.findall(r"(?<![A-Za-z])([DL]-(?:d|me|6d)?[A-Z][a-z]+[A-Za-z]*|Hex|Pen|dHex|Oct|Hept|HexA)(?![a-z])", s)

    def extractModTags(mods: str) -> List[str]:
        if pd.isna(mods) or not mods: return []
        return re.findall(r"\*([A-Za-z\d\-]+)", str(mods))


# =========================================================================
# F1: 糖链长度 × 生物活性 (Chain Length × Bioactivity Heatmap)
#
# 设计意图: 一个核心的药物化学问题 — 糖链越长,生物活性是否越强/越弱?
# 通过将 bioactivity_summary 文本分类后与 Total_Sugar_Count 交叉,
# 揭示糖链"长度-活性"相关性。
#
# Design: Core medicinal chemistry question — does sugar chain length
# correlate with bioactivity? Cross-tabulates classified bioactivity
# text with sugar count.
# =========================================================================

def _classifyBioactivity(text: str) -> List[str]:
    """Classify bioactivity summary text into standard categories.

    将自由文本生物活性摘要分类为标准类别 (抗肿瘤/抗菌/抗炎/降糖等)。
    Uses keyword matching against known bioactivity vocabulary.
    """
    if pd.isna(text) or not text:
        return []
    text = str(text).lower()
    categories = []
    keywords = {
        "Anticancer": ["anticancer", "antitumor", "cytotoxic", "antiproliferative", "tumor", "cancer"],
        "Antimicrobial": ["antimicrobial", "antibacterial", "antifungal", "antiviral", "antiparasitic"],
        "Anti-inflammatory": ["anti-inflammatory", "inflammation", "cox-2", "tnf", "nf-kb"],
        "Antidiabetic": ["antidiabetic", "hypoglycemic", "glucose", "insulin", "diabetes"],
        "Cardiovascular": ["cardiovascular", "cardioprotective", "hypotensive", "cholesterol", "lipid"],
        "Hepatoprotective": ["hepatoprotective", "liver", "hepatic"],
        "Immunomodulatory": ["immunomodulatory", "immune", "immunostimulant"],
        "Antioxidant": ["antioxidant", "radical", "oxidative"],
        "Neuroprotective": ["neuroprotective", "neurotoxic", "acetylcholinesterase"],
    }
    for category, kws in keywords.items():
        if any(kw in text for kw in kws):
            categories.append(category)
    return categories if categories else ["Other"]


def chartF1ChainLengthBioactivity(df: pd.DataFrame):
    """Heatmap: sugar chain length × bioactivity category.

    揭示"糖链越长 → 哪些生物活性更显著"的趋势。
    Reveals trends: do longer sugar chains associate with specific bioactivities?
    """
    print("[F1] Chain Length × Bioactivity Heatmap...")

    rows = []
    for _, r in df.iterrows():
        bio = r.get("bioactivity_summary", "")
        sugarCount = r.get("Total_Sugar_Count", 0)
        if pd.isna(sugarCount) or sugarCount == 0:
            continue
        categories = _classifyBioactivity(bio)
        for cat in categories:
            rows.append({"Sugar_Count": int(sugarCount), "Bioactivity": cat})

    if not rows:
        print("  [SKIP] No bioactivity data")
        return

    bdf = pd.DataFrame(rows)
    # 按 sugar count 和 bioactivity 分组
    bdf["Sugar_Bin"] = bdf["Sugar_Count"].clip(upper=8).astype(str)
    bdf.loc[bdf["Sugar_Count"] >= 8, "Sugar_Bin"] = "8+"

    topBio = bdf["Bioactivity"].value_counts().head(8).index.tolist()
    bdf = bdf[bdf["Bioactivity"].isin(topBio)]

    pivot = bdf.groupby(["Bioactivity", "Sugar_Bin"]).size().unstack(fill_value=0)
    # 排列列顺序 (Order columns logically)
    colOrder = [str(i) for i in range(1, 8)] + ["8+"]
    pivot = pivot.reindex(columns=[c for c in colOrder if c in pivot.columns], fill_value=0)

    fig = px.imshow(pivot.values, x=pivot.columns.tolist(), y=pivot.index.tolist(),
                    color_continuous_scale=HEATMAP_COLORSCALE, aspect="auto",
                    labels=dict(color="Count"), text_auto=True)
    fig.update_layout(
        title=f"Sugar Chain Length × Bioactivity (N={len(bdf):,} annotations)<br>"
              "<sub>Do longer sugar chains associate with specific bioactivities?</sub>",
        xaxis_title="Sugar Count", yaxis_title="Bioactivity Category",
    )
    saveHtml(fig, "F1_chain_vs_bioactivity.html")


# =========================================================================
# F2: 修饰类型 × ChEMBL 活性 (Modification × pChEMBL Box Plot)
#
# 设计意图: 糖修饰 (如乙酰化/硫酸化) 是否与更强的生物活性 (更高的
# pChEMBL 值) 相关? 箱线图直观展示各修饰类型的活性分布。
#
# Design: Do sugar modifications (acetylation, sulfation) correlate with
# stronger bioactivity (higher pChEMBL)? Box plots show activity distribution
# by modification type.
# =========================================================================

def chartF2ModVsPchembl(df: pd.DataFrame):
    """Box plot: modification type vs best pChEMBL value.

    展示不同糖修饰是否与更高/更低的生物活性定量指标相关。
    Shows whether specific sugar modifications correlate with quantitative bioactivity.
    """
    print("[F2] Modification × pChEMBL Box Plot...")

    subDf = df[df["ChEMBL_Best_pChEMBL"].notna()].copy()
    subDf["ChEMBL_Best_pChEMBL"] = pd.to_numeric(subDf["ChEMBL_Best_pChEMBL"], errors="coerce")
    subDf = subDf[subDf["ChEMBL_Best_pChEMBL"].notna()]

    if subDf.empty or len(subDf) < 10:
        print("  [SKIP] Insufficient pChEMBL data")
        return

    rows = []
    for _, r in subDf.iterrows():
        pchembl = r["ChEMBL_Best_pChEMBL"]
        mods = r.get("Glycan_Modifications", "")
        tags = extractModTags(str(mods)) if pd.notna(mods) else []
        if not tags:
            rows.append({"Modification": "None", "pChEMBL": pchembl})
        else:
            for tag in set(tags):
                rows.append({"Modification": tag, "pChEMBL": pchembl})

    mdf = pd.DataFrame(rows)
    topMods = mdf["Modification"].value_counts().head(10).index.tolist()
    mdf = mdf[mdf["Modification"].isin(topMods)]

    fig = px.box(mdf, x="Modification", y="pChEMBL", color="Modification",
                 color_discrete_sequence=QUAL_COLORS)
    fig.update_layout(
        title=f"Sugar Modification × pChEMBL Value (N={len(mdf):,})<br>"
              "<sub>Higher pChEMBL = stronger measured bioactivity | Box = quartiles</sub>",
        xaxis_title="Modification Type", yaxis_title="Best pChEMBL Value",
        showlegend=False,
    )
    saveHtml(fig, "F2_mod_vs_pchembl.html")


# =========================================================================
# F3: 糖类偏好 × 靶点蛋白 (Sugar Type × Target Protein Heatmap)
#
# 设计意图: 特定单糖出现时,对应哪些靶蛋白更多? 这可以揭示
# 糖结构与蛋白识别之间的偏好模式。
#
# Design: When specific monosaccharides are present, which drug targets
# are more commonly reported? Reveals sugar-protein recognition preferences.
# =========================================================================

def chartF3SugarVsTarget(df: pd.DataFrame):
    """Heatmap: top monosaccharides × top ChEMBL targets.

    特定单糖是否与特定药物靶点有统计学上的共现偏好?
    Do specific monosaccharides co-occur preferentially with certain drug targets?
    """
    print("[F3] Sugar Type × Target Protein Heatmap...")

    rows = []
    for _, r in df.iterrows():
        targets = r.get("ChEMBL_Targets", "")
        if pd.isna(targets) or not targets or str(targets) == "nan":
            continue
        sugars = r.get("Mono_List", [])
        if not isinstance(sugars, list):
            sugars = extractMonosaccharideList(str(r.get("Sugar_Sequence", "")))
        if not sugars:
            continue

        # 解析靶点列表 (Parse target list — typically pipe/comma separated)
        targetList = re.split(r"[|,;]", str(targets))
        for target in targetList:
            target = target.strip()[:40]
            if not target:
                continue
            for sugar in set(sugars):
                rows.append({"Sugar": sugar, "Target": target})

    if not rows:
        print("  [SKIP] No sugar-target co-occurrence data")
        return

    tdf = pd.DataFrame(rows)
    topSugars = tdf["Sugar"].value_counts().head(10).index.tolist()
    topTargets = tdf["Target"].value_counts().head(10).index.tolist()
    tdf = tdf[tdf["Sugar"].isin(topSugars) & tdf["Target"].isin(topTargets)]

    pivot = tdf.groupby(["Sugar", "Target"]).size().unstack(fill_value=0)
    pivot = pivot.reindex(index=topSugars, fill_value=0)

    fig = px.imshow(pivot.values, x=pivot.columns.tolist(), y=pivot.index.tolist(),
                    color_continuous_scale=HEATMAP_COLORSCALE, aspect="auto",
                    labels=dict(color="Count"), text_auto=True)
    fig.update_layout(
        title=f"Monosaccharide × Drug Target (N={len(tdf):,} co-occurrences)<br>"
              "<sub>Do specific sugars preferentially appear in molecules targeting specific proteins?</sub>",
        xaxis_title="ChEMBL Target", yaxis_title="Monosaccharide",
        margin=dict(b=150),
    )
    fig.update_xaxes(tickangle=45)
    saveHtml(fig, "F3_sugar_vs_target.html")


# =========================================================================
# F4: Lipinski RO5 违规 × 糖链长度 (RO5 Violations vs Sugar Count)
#
# 设计意图: 糖链越长 → 分子越大 → HBD/HBA 越多 → RO5 违规越多?
# 量化糖化对类药性的"稀释效应"。
#
# Design: Longer sugar chains → larger molecules → more HBD/HBA →
# more RO5 violations? Quantifies glycosylation's dilution of drug-likeness.
# =========================================================================

def chartF4Ro5VsSugarCount(df: pd.DataFrame):
    """Stacked bar: Lipinski RO5 violation distribution stratified by sugar count.

    量化糖化程度与 Lipinski 五规则合规性的关系。
    Quantifies relationship between glycosylation degree and Lipinski RO5 compliance.
    """
    print("[F4] Lipinski RO5 Violations × Sugar Count...")

    subDf = df[["Total_Sugar_Count", "lipinski_rule_of_five_violations"]].dropna().copy()
    subDf["Total_Sugar_Count"] = subDf["Total_Sugar_Count"].astype(int)
    subDf["lipinski_rule_of_five_violations"] = pd.to_numeric(
        subDf["lipinski_rule_of_five_violations"], errors="coerce"
    ).fillna(0).astype(int)
    subDf["Sugar_Bin"] = subDf["Total_Sugar_Count"].clip(upper=7)
    subDf.loc[subDf["Total_Sugar_Count"] >= 7, "Sugar_Bin"] = 7
    subDf["Violations"] = subDf["lipinski_rule_of_five_violations"].clip(upper=4).astype(str)
    subDf.loc[subDf["lipinski_rule_of_five_violations"] >= 4, "Violations"] = "4+"

    pivot = subDf.groupby(["Sugar_Bin", "Violations"]).size().unstack(fill_value=0)
    violOrder = ["0", "1", "2", "3", "4+"]
    pivot = pivot.reindex(columns=[v for v in violOrder if v in pivot.columns], fill_value=0)
    # 列归一化 (Column normalize)
    pivotPct = pivot.div(pivot.sum(axis=1), axis=0) * 100

    violColors = {"0": "#4C78A8", "1": "#72B7B2", "2": "#F2C057",
                  "3": "#E45756", "4+": "#9C1C1C"}

    fig = go.Figure()
    for viol in pivotPct.columns:
        xLabels = [f"{int(b)}" if b < 7 else "7+" for b in pivotPct.index]
        fig.add_trace(go.Bar(
            x=xLabels, y=pivotPct[viol].values,
            name=f"{viol} violation{'s' if viol != '1' else ''}",
            marker_color=violColors.get(viol, "#888"),
        ))

    fig.update_layout(
        barmode="stack",
        title=f"Lipinski RO5 Violations by Sugar Count (N={len(subDf):,})<br>"
              "<sub>100% stacked bars | More sugars → higher RO5 violation rates</sub>",
        xaxis_title="Total Sugar Count", yaxis_title="% of Molecules",
        legend_title="RO5 Violations",
        height=550,
    )
    saveHtml(fig, "F4_ro5_vs_sugar_count.html", height=550)


# =========================================================================
# F5: QED 药物相似性 × 皂苷类型 (QED × Saponin Type Violin)
#
# 设计意图: 甾体类 vs 三萜类皂苷,哪一类更接近药物? QED 综合评分
# 比单一的 RO5 更全面地反映类药性。
#
# Design: Steroidal vs Triterpenoid saponins — which class is more
# drug-like? QED score is more comprehensive than RO5 alone.
# =========================================================================

def chartF5QedBySaponinType(df: pd.DataFrame):
    """Violin plot: QED drug-likeness score by saponin sub-type.

    比较 SD 型和 TD 型皂苷的类药性分布, 指导药物先导化合物选择。
    Compares drug-likeness distributions between steroidal and triterpenoid
    saponins, guiding lead compound selection.
    """
    print("[F5] QED Drug-Likeness × Saponin Type...")

    subDf = df[df["qed_drug_likeliness"].notna()].copy()
    subDf["qed_drug_likeliness"] = pd.to_numeric(subDf["qed_drug_likeliness"], errors="coerce")
    subDf = subDf[subDf["qed_drug_likeliness"].notna()]
    subDf = subDf[subDf["Saponin_Type"].isin(["Steroidal", "Triterpenoid"])]

    if len(subDf) < 20:
        print("  [SKIP] Insufficient QED data")
        return

    fig = px.violin(subDf, x="Saponin_Type", y="qed_drug_likeliness",
                    color="Saponin_Type",
                    color_discrete_map={"Steroidal": "#4C78A8", "Triterpenoid": "#E45756"},
                    box=True, points="outliers")
    fig.update_layout(
        title=f"QED Drug-Likeness: Steroidal vs Triterpenoid (N={len(subDf):,})<br>"
              "<sub>QED integrates MW, logP, HBD/HBA, PSA, RotBonds, ArRings | Higher = more drug-like</sub>",
        xaxis_title="Saponin Type", yaxis_title="QED Score",
        showlegend=False,
    )
    saveHtml(fig, "F5_qed_by_saponin_type.html")


# =========================================================================
# Standalone execution
# =========================================================================

if __name__ == "__main__":
    print("=" * 70)
    print("  GlycoNP Saponin — F-Series: Bioactivity-Oriented Analytics")
    print("=" * 70)
    df = pd.read_csv(CSV_PATH, low_memory=False)
    totalRaw = len(df)
    for col in ["Total_Sugar_Count", "Max_Chain_Length", "exact_molecular_weight",
                "alogp", "topological_polar_surface_area", "ChEMBL_Best_pChEMBL",
                "qed_drug_likeliness"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df[df["Total_Sugar_Count"] > 0].copy()
    if "Consensus_Sugar_Sequence" in df.columns and "Sugar_Sequence" not in df.columns:
        df["Sugar_Sequence"] = df["Consensus_Sugar_Sequence"]
    df["Mono_List"] = df["Sugar_Sequence"].apply(extractMonosaccharideList)
    df["Mod_Tags"] = df["Glycan_Modifications"].apply(extractModTags)
    df["Organism_Type"] = df["Organism_Type"].fillna("Unknown")
    def classifySaponinType(row):
        sc = str(row.get("Super_Scaffold_Class", ""))
        nc = str(row.get("np_classifier_superclass", ""))
        combined = (sc + " " + nc).lower()
        if "steroid" in combined: return "Steroidal"
        elif "triterpen" in combined: return "Triterpenoid"
        return "Other"
    df["Saponin_Type"] = df.apply(classifySaponinType, axis=1)
    print(f"  Loaded {len(df):,} saponins\n")

    import traceback
    for func in [chartF1ChainLengthBioactivity, chartF2ModVsPchembl,
                 chartF3SugarVsTarget, chartF4Ro5VsSugarCount, chartF5QedBySaponinType]:
        try:
            func(df)
        except Exception as e:
            print(f"  [ERROR] {func.__name__}: {e}")
            traceback.print_exc()
    print("\n  F-Series complete!")
