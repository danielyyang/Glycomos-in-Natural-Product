#!/usr/bin/env python3
"""
==========================================================================
  [EN] GlycoNP Saponin Database — Synthesis-Oriented Analytics (E-Series)
       Charts designed to guide synthetic chemists: frequent building blocks,
       aglycone glycosylation hotspots, inter-glycan linkage distribution,
       branching degree analysis, and 1,2-cis glycoside prevalence.

  [CN] GlycoNP 皂苷子库 — 合成导向分析图谱 (E 系列)
       为合成化学家设计的专题图表：高频合成砌块、苷元糖基化热点位、
       糖链内部连接分布、分支度分析、1,2-顺式糖苷键占比。

  Input:  reports/GlycoNP_Saponin_DB.csv
  Output: reports/saponin_figures/E*.html + E*.png
  [TEST DATA ONLY]
==========================================================================
"""
import os
import re
import sys
import json
import warnings
from pathlib import Path
from collections import Counter
from typing import List, Dict, Tuple

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

warnings.filterwarnings("ignore")

# ── 全局变量: 由 _loadModule 注入或独立运行时自行定义 ─────────────────
# Globals: injected by _loadModule, or self-defined in standalone mode
try:
    _ = LAYOUT_DEFAULTS  # type: ignore
except NameError:
    # 独立运行模式 (Standalone mode)
    BASE_DIR = Path(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
    CSV_PATH = BASE_DIR / "reports" / "GlycoNP_Saponin_DB.csv"
    OUT_DIR = BASE_DIR / "reports" / "saponin_figures"
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    LAYOUT_DEFAULTS = dict(
        template="plotly_white",
        font=dict(family="Arial", size=13, color="#222"),
        paper_bgcolor="white",
        plot_bgcolor="white",
        margin=dict(l=80, r=40, t=80, b=80),
        title_font=dict(size=17, family="Arial", color="#222"),
    )
    HEATMAP_COLORSCALE = "Viridis"
    QUAL_COLORS = px.colors.qualitative.Set2 + px.colors.qualitative.Pastel1

    def saveHtml(fig, filename: str, height: int = 650, width: int = 1100):
        """Save Plotly figure as HTML + PNG (standalone fallback)."""
        curMargin = fig.layout.margin
        customL = getattr(curMargin, "l", None)
        customR = getattr(curMargin, "r", None)
        customT = getattr(curMargin, "t", None)
        customB = getattr(curMargin, "b", None)
        fig.update_layout(**LAYOUT_DEFAULTS)
        fig.update_layout(margin=dict(
            l=customL if customL is not None else 80,
            r=customR if customR is not None else 40,
            t=customT if customT is not None else 80,
            b=customB if customB is not None else 80,
        ))
        if height:
            fig.update_layout(height=height)
        if width:
            fig.update_layout(width=width)
        outPath = OUT_DIR / filename
        fig.write_html(str(outPath), include_plotlyjs="cdn")
        print(f"  -> Saved: {outPath.name}")
        pngName = filename.replace(".html", ".png")
        pngPath = OUT_DIR / pngName
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
            print(f"  -> Saved: {pngName}")
        except Exception as e:
            print(f"  [WARN] PNG export failed ({e}).")

    def extractMonosaccharideList(seq: str) -> List[str]:
        """Extract monosaccharide names from a Sugar_Sequence string."""
        if pd.isna(seq) or not seq:
            return []
        s = re.sub(r"\([ab]\d-\d\)", " ", str(seq))
        s = s.replace(";", " ").replace("[", " ").replace("]", " ")
        pattern = r"(?<![A-Za-z])([DL]-(?:d|me|6d)?[A-Z][a-z]+[A-Za-z]*|Hex|Pen|dHex|Oct|Hept|HexA)(?![a-z])"
        return re.findall(pattern, s)

    def extractModTags(mods: str) -> List[str]:
        """Extract modification tags from Glycan_Modifications."""
        if pd.isna(mods) or not mods:
            return []
        return re.findall(r"\*([A-Za-z\d\-]+)", str(mods))

    def parseSequenceToEdges(seq: str) -> List[Dict]:
        """Parse Sugar_Sequence into structured edges."""
        if pd.isna(seq) or not seq:
            return []
        edges = []
        pattern = r"([A-Za-z\d\-]+)-\(([ab])(\d)-(\d)\)-([A-Za-z\d\-]+)"
        for match in re.finditer(pattern, str(seq)):
            donor, anomerChar, posFrom, posTo, acceptor = match.groups()
            anomer = "α" if anomerChar == "a" else "β"
            edges.append({"donor": donor, "acceptor": acceptor,
                          "anomer": anomer, "pos": f"{posFrom}→{posTo}"})
        return edges


# =========================================================================
# E1: 连接感知二糖砌块频率 (Linkage-Aware Disaccharide Synthon Frequency)
#
# 设计意图: S1 只做了"文本滑窗"统计，未区分连接方式。
# 对合成化学家而言，Rha(α1→2)Glc 和 Rha(β1→4)Glc 是完全不同的砌块，
# 需要不同的合成方法。本图以 "Donor(anomer pos→pos)Acceptor" 四元组
# 为唯一键统计频率，直接回答"实验室应预制哪些二糖砌块"的问题。
#
# Design: S1 uses text sliding-window without differentiating linkage types.
# For synthetic chemists, Rha(α1→2)Glc vs Rha(β1→4)Glc require entirely
# different synthetic methods. This chart uses Donor(anomer-pos)→Acceptor
# as the unique key, directly answering "which disaccharide building blocks
# should be pre-synthesized in the lab?"
# =========================================================================

def chartE1LinkageAwareSynthons(df: pd.DataFrame):
    """Bar chart of top 30 linkage-aware disaccharide synthons.

    解析每条糖链序列中的二糖连接四元组 (Donor-anomer-position-Acceptor)
    并统计全库频率, 为合成实验室的"通用砌块"预制提供数据依据。
    Parses disaccharide linkage quadruples from every sugar sequence and
    counts global frequency, guiding pre-synthesis of universal building blocks.
    """
    print("[E1] Linkage-Aware Disaccharide Synthon Frequency...")
    allEdges: List[str] = []
    for seq in df["Sugar_Sequence"].dropna():
        edges = parseSequenceToEdges(str(seq))
        for e in edges:
            # 构建四元组标签 (Build quadruple label)
            # 格式: Donor-(anomer pos→pos)-Acceptor
            # 设计意图: 完整保留连接信息, 不丢失任何合成相关细节
            label = f"{e['donor']}-({e['anomer']}{e['pos']})-{e['acceptor']}"
            allEdges.append(label)

    if not allEdges:
        print("  [SKIP] No linkage edges found")
        return

    edgeCounts = Counter(allEdges).most_common(30)
    labels, counts = zip(*edgeCounts)

    # 颜色编码: α 用暖色, β 用冷色 (Color: warm=α, cool=β)
    # 设计意图: 一眼区分异头碳构型, 方便合成化学家判断合成难度
    colors = []
    for lbl in labels:
        if "α" in lbl:
            colors.append("#E45756")  # 暖色 → α
        elif "β" in lbl:
            colors.append("#4C78A8")  # 冷色 → β
        else:
            colors.append("#72B7B2")

    fig = go.Figure(go.Bar(
        y=list(reversed(labels)), x=list(reversed(counts)),
        orientation="h",
        marker_color=list(reversed(colors)),
        text=list(reversed(counts)), textposition="outside",
    ))
    fig.update_layout(
        title=f"Top 30 Linkage-Aware Disaccharide Synthons (N={sum(counts):,})<br>"
              "<sub>Red=α-linked, Blue=β-linked | Each bar = one unique building block for synthesis</sub>",
        xaxis_title="Occurrence Count",
        yaxis_title="Disaccharide Synthon (Donor-Linkage-Acceptor)",
        height=900,
        margin=dict(l=250),
    )
    saveHtml(fig, "E1_linkage_aware_synthons.html", height=900)


# =========================================================================
# E2: 苷元碳位修饰热图 (Aglycone Carbon Position × Sugar Type Heatmap)
#
# 设计意图: 合成皂苷时, 需要知道骨架上哪些碳位最常被糖基化。
# 例如齐墩果烷型皂苷 C-3 多连中性糖、C-28 多连糖醛酸,
# 据此可设计正交保护基组合。本图从 Glycan-Aglycone_Bond_Detail JSON
# 中解析碳位信息, 生成 Position × Sugar 的交叉热力图。
#
# Design: When synthesizing saponins, chemists need to know which aglycone
# carbon positions are most frequently glycosylated. This chart parses
# Glycan-Aglycone_Bond_Detail JSON for position info and creates a
# Position × Sugar heatmap, guiding orthogonal protection group design.
# =========================================================================

def chartE2AglyconePositionHeatmap(df: pd.DataFrame):
    """Heatmap: which sugars are attached at which aglycone carbon positions.

    从 Bond_Detail JSON 中提取苷元碳位编号 (如 C-3, C-28),
    交叉糖类型生成热力图。直接指导保护基策略设计。
    Extracts aglycone carbon positions from Bond_Detail JSON,
    cross-tabulates with sugar type for a heatmap guiding protection strategy.
    """
    print("[E2] Aglycone Carbon Position × Sugar Type Heatmap...")

    records: List[Dict] = []
    for _, row in df.iterrows():
        bondRaw = row.get("Glycan-Aglycone_Bond_Detail", "[]")
        try:
            bonds = json.loads(str(bondRaw)) if pd.notna(bondRaw) and str(bondRaw) not in ("nan", "", "[]") else []
        except Exception:
            bonds = []

        for bd in bonds:
            sugar = bd.get("sugar", "?")
            target = bd.get("target", "")
            bond = bd.get("bond", "")
            if "Aglycon" not in target and "aglycon" not in target.lower():
                continue

            # 从 bond 字段提取连接类型 (Extract linkage type from bond field)
            # 格式通常为 "α-O-linked", "β-O-linked" 等
            linkElement = "O"
            if "-N-" in bond:
                linkElement = "N"
            elif "-C-" in bond:
                linkElement = "C"
            elif "-S-" in bond:
                linkElement = "S"

            records.append({
                "Sugar": sugar,
                "Linkage_Element": linkElement,
                "Saponin_Type": row.get("Saponin_Type", "Unknown"),
            })

    if not records:
        print("  [SKIP] No aglycone bond data")
        return

    rdf = pd.DataFrame(records)

    # -- Panel A: 连接根糖 × 连接元素 (Root Sugar × Linkage Element) --
    topSugars = rdf["Sugar"].value_counts().head(15).index.tolist()
    sub = rdf[rdf["Sugar"].isin(topSugars)]
    pivot = sub.groupby(["Sugar", "Linkage_Element"]).size().unstack(fill_value=0)
    pivot = pivot.reindex(index=topSugars, fill_value=0)

    fig = go.Figure(data=go.Heatmap(
        z=pivot.values,
        x=pivot.columns.tolist(), y=pivot.index.tolist(),
        colorscale=HEATMAP_COLORSCALE,
        text=[[str(int(v)) for v in row] for row in pivot.values],
        texttemplate="%{text}", textfont=dict(size=12),
        colorbar=dict(title="Count"),
    ))
    fig.update_layout(
        title=f"Root Sugar × Aglycone Linkage Element (N={len(rdf):,})<br>"
              "<sub>Which sugars connect to the aglycone via O/N/C/S bonds?</sub>",
        xaxis_title="Linkage Element", yaxis_title="Root Sugar",
        height=550,
    )
    saveHtml(fig, "E2a_root_sugar_linkage.html", height=550)

    # -- Panel B: 根糖 × 皂苷类型 (Root Sugar × Saponin Type) --
    # 设计意图: SD 型 vs TD 型皂苷的糖基化偏好差异
    # Design: Steroidal vs Triterpenoid glycosylation preference comparison
    sub2 = rdf[rdf["Sugar"].isin(topSugars) & rdf["Saponin_Type"].isin(["Steroidal", "Triterpenoid"])]
    if not sub2.empty:
        pivot2 = sub2.groupby(["Sugar", "Saponin_Type"]).size().unstack(fill_value=0)
        colTotals = pivot2.sum(axis=0)
        pivotPct = pivot2.div(colTotals, axis=1) * 100
        pivotPct.columns = [f"{c} (N={int(colTotals[c]):,})" for c in pivotPct.columns]
        pivotPct = pivotPct.reindex(index=topSugars, fill_value=0)

        fig2 = px.imshow(pivotPct.values, x=pivotPct.columns.tolist(),
                         y=pivotPct.index.tolist(),
                         color_continuous_scale=HEATMAP_COLORSCALE, aspect="auto",
                         labels=dict(color="Col %"), text_auto=".1f")
        fig2.update_layout(
            title="Root Sugar × Saponin Type (Column-Normalized)<br>"
                  "<sub>Steroidal vs Triterpenoid: which sugars anchor to the aglycone?</sub>",
            xaxis_title="Saponin Type", yaxis_title="Root Sugar",
        )
        saveHtml(fig2, "E2b_root_sugar_saponin_type.html")


# =========================================================================
# E3: 糖链内部连接矩阵 (Inter-Glycan Linkage Distribution Matrix)
#
# 设计意图: A8 系列统计了糖-苷元连接键, 但缺少糖-糖之间的内部连接统计。
# 合成化学家需要知道 α(1→2), β(1→3), β(1→4), α(1→6) 等内部键的精确比例,
# 因为不同键型需要完全不同的合成方法学 (如 1,2-顺式键需特殊催化剂)。
#
# Design: A8 series covers sugar-aglycone bonds, but lacks intra-glycan linkage
# statistics. Synthetically, α(1→2) vs β(1→4) require entirely different
# methodologies. This chart shows the precise distribution of all internal
# glycan-glycan linkage types.
# =========================================================================

def chartE3InterGlycanLinkageMatrix(df: pd.DataFrame):
    """Heatmap: anomer (α/β) × linkage position (1→2, 1→3, etc.) for intra-glycan bonds.

    从糖链序列中提取所有内部糖-糖连接的异头碳构型和连接位置,
    生成 α/β × Position 的交叉频率矩阵。
    Extracts all inter-glycan linkage info from sugar sequences,
    creating an anomer × position frequency matrix.
    """
    print("[E3] Inter-Glycan Linkage Distribution Matrix...")

    allEdges = []
    for seq in df["Sugar_Sequence"].dropna():
        edges = parseSequenceToEdges(str(seq))
        allEdges.extend(edges)

    if not allEdges:
        print("  [SKIP] No edges found")
        return

    edf = pd.DataFrame(allEdges)

    # -- Panel A: α/β × Position 热力图 (Heatmap) --
    pivot = edf.groupby(["anomer", "pos"]).size().unstack(fill_value=0)
    # 按位置排序 (Sort by position label)
    sortedCols = sorted(pivot.columns.tolist(), key=lambda x: (int(x[0]), int(x[-1])))
    pivot = pivot.reindex(columns=sortedCols, fill_value=0)

    fig = go.Figure(data=go.Heatmap(
        z=pivot.values,
        x=pivot.columns.tolist(), y=pivot.index.tolist(),
        colorscale=HEATMAP_COLORSCALE,
        text=[[f"{int(v):,}" for v in row] for row in pivot.values],
        texttemplate="%{text}", textfont=dict(size=13),
        colorbar=dict(title="Count"),
    ))
    fig.update_layout(
        title=f"Inter-Glycan Linkage Matrix: Anomer × Position (N={len(edf):,} bonds)<br>"
              "<sub>All sugar↔sugar bonds in 26K saponins | Higher counts = more common linkage</sub>",
        xaxis_title="Linkage Position", yaxis_title="Anomeric Configuration",
        height=400,
    )
    saveHtml(fig, "E3a_interglycan_linkage_matrix.html", height=400)

    # -- Panel B: Donor-Acceptor 组合频率 (Top pairs) --
    # 设计意图: 哪些具体的"供体→受体"糖对最常见？
    # Design: Which specific donor→acceptor sugar pairs are most common?
    edf["Pair"] = edf["donor"] + "→" + edf["acceptor"]
    pairCounts = edf["Pair"].value_counts().head(20)

    fig2 = go.Figure(go.Bar(
        y=pairCounts.index[::-1], x=pairCounts.values[::-1],
        orientation="h", marker_color=QUAL_COLORS[:len(pairCounts)],
        text=pairCounts.values[::-1], textposition="outside",
    ))
    fig2.update_layout(
        title=f"Top 20 Sugar→Sugar Pairs (N={pairCounts.sum():,})<br>"
              "<sub>Most frequent donor→acceptor combinations across all glycan chains</sub>",
        xaxis_title="Occurrence Count",
        yaxis_title="Donor → Acceptor",
        height=650,
        margin=dict(l=200),
    )
    saveHtml(fig2, "E3b_sugar_pair_frequency.html", height=650)


# =========================================================================
# E4: 单糖分支度分布 (Monosaccharide Branching Degree Distribution)
#
# 设计意图: 一个单糖上同时被多少个其他糖基化 (如 3,4-位同时被修饰)
# 的频率直接决定了合成难度。空间位阻极大的多分支糖链在合成中
# 极易失败, 本图帮助评估"合成死胡同"的风险。
#
# Design: How many glycan chains branch off a single monosaccharide?
# Dual/triple-glycosylation at adjacent positions (e.g. 3,4-di-glycosylation)
# creates high steric hindrance and often leads to synthetic failure.
# This chart maps the "synthetic dead-end" risk.
# =========================================================================

def chartE4BranchingDegree(df: pd.DataFrame):
    """Distribution of monosaccharide branching degree (how many sugars branch
    from a single residue).

    通过解析分支序列中的 [...] 语法，统计每个分支点挂载的子链数。
    Parses branch notation [...] to count sub-chains per branch point.
    """
    print("[E4] Monosaccharide Branching Degree Distribution...")

    branchDegrees: List[int] = []
    branchSugars: List[str] = []

    for seq in df["Sugar_Sequence"].dropna():
        seqStr = str(seq)
        if "[" not in seqStr:
            continue

        # 解析分支结构 (Parse branch structure)
        # 每个 [...] 代表一个分支; 连续多个 [...] 意味着同一个糖上的多分支
        # Each [...] = one branch; consecutive [...] = multi-branched sugar
        chains = seqStr.split(";")
        for chain in chains:
            chain = chain.strip()
            if "[" not in chain:
                continue

            # 统计分支点数量 (Count branch points)
            # 简化策略: 通过 ] 后紧跟的糖名来识别分支点
            # Simplified: after each ], the next sugar name is the branch point
            # 分支度 = 该糖直接连接的 [...] 块数 + 1 (主链方向)
            branchBlocks = re.findall(r"\[([^\]]+)\]", chain)
            if branchBlocks:
                # 分支度 = 分支数 + 1 (主链方向本身也算一个连接)
                # Branching degree = number of branches + 1 (main chain)
                degree = len(branchBlocks) + 1
                branchDegrees.append(degree)

                # 尝试识别分支点糖名 (Identify branch-point sugar name)
                # 模式: ]-SugarName 后面的第一个糖
                afterBracket = re.findall(r"\]-([A-Za-z\d\-]+)", chain)
                for sugar in afterBracket:
                    branchSugars.append(sugar)

    if not branchDegrees:
        print("  [SKIP] No branch data")
        return

    # -- Panel A: 分支度分布直方图 (Branching degree histogram) --
    degreeCounts = Counter(branchDegrees)
    degrees = sorted(degreeCounts.keys())
    counts = [degreeCounts[d] for d in degrees]

    fig = go.Figure(go.Bar(
        x=[f"Degree {d}" for d in degrees], y=counts,
        marker_color=QUAL_COLORS[:len(degrees)],
        text=counts, textposition="outside",
    ))
    fig.update_layout(
        title=f"Branching Degree Distribution (N={sum(counts):,} branch points)<br>"
              "<sub>Degree 2 = one branch + main chain; Degree 3+ = multi-branch (high steric hindrance)</sub>",
        xaxis_title="Branching Degree",
        yaxis_title="Count (Branch Points)",
        height=500,
    )
    saveHtml(fig, "E4a_branching_degree_dist.html", height=500)

    # -- Panel B: 分支点糖类型分布 (Branch point sugar types) --
    if branchSugars:
        sugarCounts = Counter(branchSugars).most_common(15)
        names, cnts = zip(*sugarCounts)
        fig2 = go.Figure(go.Bar(
            y=list(reversed(names)), x=list(reversed(cnts)),
            orientation="h", marker_color="#B279A2",
            text=list(reversed(cnts)), textposition="outside",
        ))
        fig2.update_layout(
            title=f"Most Common Branch-Point Sugars (N={sum(cnts):,})<br>"
                  "<sub>These sugars most often serve as multi-branching hubs in glycan chains</sub>",
            xaxis_title="Count", yaxis_title="Branch-Point Sugar",
            height=550,
            margin=dict(l=180),
        )
        saveHtml(fig2, "E4b_branch_point_sugars.html", height=550)


# =========================================================================
# E5: 1,2-顺式糖苷键占比分析 (1,2-cis Glycoside Prevalence)
#
# 设计意图: 1,2-顺式糖苷键 (α-Glc, β-Man, α-GalNAc 等) 是合成化学
# 最大的挑战之一, 需要手性辅助基团或特殊催化剂。如果数据库中
# 大量存在此类键, 则项目需要额外的合成方法学储备。
#
# Design: 1,2-cis glycosidic bonds (α-Glc, β-Man, α-GalNAc) are among
# the hardest synthetic challenges, requiring chiral auxiliaries or
# special catalysts. If prevalent in the database, the synthesis project
# needs additional methodology investment.
# =========================================================================

# 1,2-cis 和 1,2-trans 的定义参考 Demchenko 2008:
# D-Gluco/D-Galacto 系列: α = 1,2-cis, β = 1,2-trans
# D-Manno 系列: β = 1,2-cis, α = 1,2-trans
# Reference: Demchenko 2008 for 1,2-cis/trans definitions.
CIS_12_RULES: Dict[str, str] = {
    "D-Glc":    "α",  # α-D-Glc = 1,2-cis | β-D-Glc = 1,2-trans (自然界最常见)
    "D-Gal":    "α",  # α-D-Gal = 1,2-cis
    "D-GlcNAc": "α",  # α-D-GlcNAc = 1,2-cis
    "D-GalNAc": "α",  # α-D-GalNAc = 1,2-cis
    "D-GlcA":   "α",  # α-D-GlcA = 1,2-cis
    "D-GalA":   "α",  # α-D-GalA = 1,2-cis
    "D-Man":    "β",  # β-D-Man = 1,2-cis | α-D-Man = 1,2-trans
    "D-ManNAc": "β",  # β-D-ManNAc = 1,2-cis
    "L-Rha":    "β",  # β-L-Rha = 1,2-cis (6-deoxy-L-Man 系列)
    "L-Fuc":    "β",  # β-L-Fuc = 1,2-cis (6-deoxy-L-Gal 系列)
}

def chartE5CisTransPrevalence(df: pd.DataFrame):
    """Stacked bar + pie: proportion of 1,2-cis vs 1,2-trans glycosidic bonds.

    基于 Demchenko 2008 定义的 1,2-cis/trans 规则, 对全库所有内部
    糖苷键进行分类, 量化合成化学的"难度负荷"。
    Using Demchenko 2008 1,2-cis/trans rules, classifies all intra-glycan
    bonds to quantify the "synthesis difficulty load".
    """
    print("[E5] 1,2-cis / 1,2-trans Glycoside Prevalence...")

    cisCount = 0
    transCount = 0
    unknownCount = 0
    cisByDonor: Dict[str, int] = Counter()
    transByDonor: Dict[str, int] = Counter()

    for seq in df["Sugar_Sequence"].dropna():
        edges = parseSequenceToEdges(str(seq))
        for e in edges:
            donor = e["donor"]
            anomer = e["anomer"]

            # 查找 1,2-cis 规则 (Look up 1,2-cis rule for this donor)
            cisAnomer = CIS_12_RULES.get(donor)
            if cisAnomer is None:
                unknownCount += 1
                continue

            if anomer == cisAnomer:
                cisCount += 1
                cisByDonor[donor] += 1
            else:
                transCount += 1
                transByDonor[donor] += 1

    total = cisCount + transCount + unknownCount
    if total == 0:
        print("  [SKIP] No linkage data")
        return

    # -- Panel A: 总体 Pie (Overall pie chart) --
    fig = make_subplots(
        rows=1, cols=2,
        specs=[[{"type": "pie"}, {"type": "bar"}]],
        subplot_titles=["Overall 1,2-cis vs 1,2-trans", "1,2-cis Bonds by Donor Sugar"],
    )

    fig.add_trace(go.Pie(
        labels=["1,2-cis (Hard)", "1,2-trans (Easy)", "Unknown"],
        values=[cisCount, transCount, unknownCount],
        marker_colors=["#E45756", "#4C78A8", "#CCCCCC"],
        textinfo="label+percent",
        hole=0.35,
    ), row=1, col=1)

    # -- Panel B: 按供体糖分类的 1,2-cis 键 (1,2-cis bonds by donor) --
    cisSorted = sorted(cisByDonor.items(), key=lambda x: -x[1])
    if cisSorted:
        names, cnts = zip(*cisSorted[:12])
        fig.add_trace(go.Bar(
            x=list(names), y=list(cnts),
            marker_color="#E45756",
            text=list(cnts), textposition="outside",
            showlegend=False,
        ), row=1, col=2)

    fig.update_layout(
        title=f"1,2-cis vs 1,2-trans Glycoside Distribution (N={total:,} bonds)<br>"
              "<sub>1,2-cis bonds (red) require chiral auxiliaries or special catalysts for synthesis</sub>",
        height=500,
    )
    saveHtml(fig, "E5_cis_trans_prevalence.html", height=500)


# =========================================================================
# Standalone execution entry point
# 独立运行入口: 直接生成 E 系列所有图表
# =========================================================================

if __name__ == "__main__":
    print("=" * 70)
    print("  GlycoNP Saponin — E-Series: Synthesis-Oriented Analytics")
    print("=" * 70)

    df = pd.read_csv(CSV_PATH, low_memory=False)
    totalRaw = len(df)
    for col in ["Total_Sugar_Count", "Max_Chain_Length"]:
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
        if "steroid" in combined:
            return "Steroidal"
        elif "triterpen" in combined:
            return "Triterpenoid"
        return "Other"
    df["Saponin_Type"] = df.apply(classifySaponinType, axis=1)
    print(f"  Loaded {len(df):,} saponins (from {totalRaw:,})\n")

    for func in [chartE1LinkageAwareSynthons, chartE2AglyconePositionHeatmap,
                 chartE3InterGlycanLinkageMatrix, chartE4BranchingDegree,
                 chartE5CisTransPrevalence]:
        try:
            func(df)
        except Exception as e:
            import traceback
            print(f"  [ERROR] {func.__name__}: {e}")
            traceback.print_exc()

    print("\n  E-Series complete!")
