"""
GlycoNP 科级网络图 — 骨架-糖链-物种 交互关系图
GlycoNP Family Network — Scaffold ↔ Sugar ↔ Organism Interactive Graph

拓扑结构 (Network Topology):
  🔵 中心节点 (大): Murcko_Scaffold 苷元骨架
  🟡 一级分支 (中): Sugar_Sequence 糖链序列
  🟢 叶子节点 (小): Organism 具体物种

每个科 (Family) 生成一张独立的交互式 HTML 网络图。
直观展示：同一个骨架上接了哪几种糖链，被哪些物种使用。

使用方法 (Usage):
  python scripts/plot_family_network.py --family Solanaceae [--input PATH]
  python scripts/plot_family_network.py --family Fabaceae --min-scaffold 3
"""
import argparse
import os
import sys
from collections import Counter, defaultdict
from typing import Dict, List, Set

import pandas as pd

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))


# =====================================================================
# 1. 属名→科 映射字典 (Genus → Family Dictionary)
# =====================================================================
# 内置常见科属映射, 覆盖 COCONUT 中高频的天然产物来源科
# 数据来源: APG IV 分类系统 [Reference Only]

GENUS_TO_FAMILY: Dict[str, str] = {
    # 茄科 Solanaceae — 甾体生物碱糖苷 (Steroidal alkaloid glycosides)
    "Solanum": "Solanaceae", "Lycianthes": "Solanaceae",
    "Capsicum": "Solanaceae", "Nicotiana": "Solanaceae",
    "Datura": "Solanaceae", "Atropa": "Solanaceae",
    "Withania": "Solanaceae", "Physalis": "Solanaceae",
    "Lycium": "Solanaceae", "Petunia": "Solanaceae",
    "Cestrum": "Solanaceae", "Hyoscyamus": "Solanaceae",

    # 豆科 Fabaceae — 三萜皂苷 (Triterpenoid saponins)
    "Glycine": "Fabaceae", "Astragalus": "Fabaceae",
    "Glycyrrhiza": "Fabaceae", "Sophora": "Fabaceae",
    "Trifolium": "Fabaceae", "Medicago": "Fabaceae",
    "Pisum": "Fabaceae", "Lupinus": "Fabaceae",
    "Acacia": "Fabaceae", "Bauhinia": "Fabaceae",
    "Dalbergia": "Fabaceae", "Pueraria": "Fabaceae",

    # 菊科 Asteraceae — 倍半萜内酯 (Sesquiterpene lactones)
    "Artemisia": "Asteraceae", "Helianthus": "Asteraceae",
    "Taraxacum": "Asteraceae", "Chrysanthemum": "Asteraceae",
    "Senecio": "Asteraceae", "Stevia": "Asteraceae",
    "Calendula": "Asteraceae", "Inula": "Asteraceae",
    "Echinacea": "Asteraceae", "Achillea": "Asteraceae",

    # 五加科 Araliaceae — 人参皂苷 (Ginsenosides)
    "Panax": "Araliaceae", "Aralia": "Araliaceae",
    "Eleutherococcus": "Araliaceae", "Hedera": "Araliaceae",
    "Dendropanax": "Araliaceae",

    # 百合科/天门冬科 Asparagaceae (旧 Liliaceae)
    "Asparagus": "Asparagaceae", "Dracaena": "Asparagaceae",
    "Polygonatum": "Asparagaceae", "Convallaria": "Asparagaceae",
    "Yucca": "Asparagaceae", "Agave": "Asparagaceae",

    # 十字花科 Brassicaceae — 硫代葡萄糖苷 (Glucosinolates)
    "Brassica": "Brassicaceae", "Arabidopsis": "Brassicaceae",
    "Raphanus": "Brassicaceae", "Sinapis": "Brassicaceae",
    "Armoracia": "Brassicaceae", "Erysimum": "Brassicaceae",
    "Isatis": "Brassicaceae", "Lepidium": "Brassicaceae",

    # 蔷薇科 Rosaceae
    "Rosa": "Rosaceae", "Prunus": "Rosaceae",
    "Malus": "Rosaceae", "Rubus": "Rosaceae",
    "Fragaria": "Rosaceae", "Potentilla": "Rosaceae",
    "Crataegus": "Rosaceae", "Pyrus": "Rosaceae",

    # 唇形科 Lamiaceae
    "Salvia": "Lamiaceae", "Mentha": "Lamiaceae",
    "Thymus": "Lamiaceae", "Ocimum": "Lamiaceae",
    "Scutellaria": "Lamiaceae", "Prunella": "Lamiaceae",
    "Leonurus": "Lamiaceae", "Perilla": "Lamiaceae",

    # 毛茛科 Ranunculaceae
    "Cimicifuga": "Ranunculaceae", "Actaea": "Ranunculaceae",
    "Aconitum": "Ranunculaceae", "Clematis": "Ranunculaceae",
    "Ranunculus": "Ranunculaceae", "Helleborus": "Ranunculaceae",
    "Paeonia": "Paeoniaceae",

    # 禾本科 Poaceae
    "Oryza": "Poaceae", "Avena": "Poaceae",
    "Triticum": "Poaceae", "Zea": "Poaceae",

    # 芸香科 Rutaceae — 香豆素 (Coumarins)
    "Citrus": "Rutaceae", "Ruta": "Rutaceae",
    "Phellodendron": "Rutaceae", "Zanthoxylum": "Rutaceae",

    # 微生物
    "Streptomyces": "Streptomycetaceae",
    "Saccharopolyspora": "Pseudonocardiaceae",
    "Bacillus": "Bacillaceae",
}


def inferFamily(organism: str) -> str:
    """
    从物种全名提取属名, 查表返回科名。
    Extract genus from species name and lookup family.
    "Solanum lycopersicum" → "Solanaceae"
    """
    if not organism or str(organism) == "nan":
        return "Unknown"
    genus = str(organism).strip().split()[0]
    return GENUS_TO_FAMILY.get(genus, "Unknown")


# =====================================================================
# 2. 数据过滤与准备 (Data Filtering & Preparation)
# =====================================================================

def filterByFamily(df: pd.DataFrame, targetFamily: str) -> pd.DataFrame:
    """
    筛选属于目标科的化合物, 展开多物种行。
    Filter compounds belonging to target family, exploding multi-organism rows.

    Args:
        df: 全量管线输出
        targetFamily: 目标科名 (如 'Solanaceae')

    Returns:
        过滤后的 DataFrame, 每行对应一个 (compound, organism) 对
    """
    rows = []
    for _, row in df.iterrows():
        orgs = str(row.get("organisms", "")).split("|")
        scaffold = str(row.get("Murcko_Scaffold", ""))
        sequence = str(row.get("Sugar_Sequence", ""))
        superclass = str(row.get("Superclass", ""))
        smiles = str(row.get("canonical_smiles", ""))
        mods = str(row.get("Glycan_Modifications", ""))

        if scaffold in ("", "nan"):
            scaffold = "Aliphatic / Unknown"
        if sequence in ("", "nan"):
            sequence = "Unknown"

        for org in orgs:
            org = org.strip()
            if not org or org == "nan":
                continue
            family = inferFamily(org)
            if family == targetFamily:
                rows.append({
                    "organism": org,
                    "family": family,
                    "scaffold": scaffold,
                    "sugar_sequence": sequence,
                    "superclass": superclass,
                    "modifications": mods if mods != "nan" else "",
                    "smiles": smiles,
                })

    return pd.DataFrame(rows)


# =====================================================================
# 3. 网络图构建 (Network Graph Construction)
# =====================================================================

def buildFamilyNetwork(
    familyDf: pd.DataFrame,
    targetFamily: str,
    outputPath: str,
    minScaffoldCount: int = 2,
):
    """
    构建 Scaffold → Sugar → Organism 三层网络图。
    Build 3-tier network: Scaffold → Sugar → Organism.

    Args:
        familyDf: filterByFamily() 的输出
        targetFamily: 科名 (用于标题)
        outputPath: HTML 输出路径
        minScaffoldCount: 最少化合物数门槛 (避免孤立节点)
    """
    from pyvis.network import Network

    if familyDf.empty:
        print(f"  [WARNING] No compounds found for {targetFamily}")
        return

    # 过滤低频骨架 (Filter rare scaffolds)
    scaffoldCounts = familyDf["scaffold"].value_counts()
    validScaffolds = scaffoldCounts[scaffoldCounts >= minScaffoldCount].index.tolist()
    plotDf = familyDf[familyDf["scaffold"].isin(validScaffolds)]

    if plotDf.empty:
        print(f"  [WARNING] No scaffolds with >= {minScaffoldCount} compounds for {targetFamily}")
        # 如果太严格就放宽 (Relax threshold)
        plotDf = familyDf
        validScaffolds = familyDf["scaffold"].unique().tolist()[:15]  # 最多 15 个
        plotDf = familyDf[familyDf["scaffold"].isin(validScaffolds)]

    # 创建网络 (Create network)
    net = Network(
        height="800px",
        width="100%",
        bgcolor="#1a1a2e",
        font_color="white",
        directed=False,
        notebook=False,
    )

    # 物理引擎设置 (Physics settings — force-directed layout)
    net.force_atlas_2based(
        gravity=-80,
        central_gravity=0.008,
        spring_length=150,
        spring_strength=0.04,
        damping=0.4,
    )

    # 收集唯一节点 (Collect unique nodes)
    scaffolds: Set[str] = set()
    sugars: Set[str] = set()
    organisms: Set[str] = set()

    # 收集连接 (Collect edges)
    scaffoldSugar: Counter = Counter()   # (scaffold, sugar) → count
    sugarOrganism: Counter = Counter()   # (sugar, organism) → count
    scaffoldOrganism: Set = set()        # (scaffold, organism) — for hover info

    for _, row in plotDf.iterrows():
        sc = row["scaffold"]
        su = row["sugar_sequence"]
        org = row["organism"]

        scaffolds.add(sc)
        sugars.add(su)
        organisms.add(org)

        scaffoldSugar[(sc, su)] += 1
        sugarOrganism[(su, org)] += 1
        scaffoldOrganism.add((sc, org))

    # 骨架 SMILES 截断显示 (Truncate scaffold SMILES for labels)
    def scaffoldLabel(smi: str) -> str:
        if len(smi) > 40:
            return smi[:37] + "..."
        return smi

    # 添加骨架节点 (Add scaffold nodes — large, blue)
    for sc in scaffolds:
        count = int(scaffoldCounts.get(sc, 0))
        size = int(max(25, min(55, 15 + count * 3)))
        label = scaffoldLabel(sc)
        # 获取该骨架的大类
        scDf = plotDf[plotDf["scaffold"] == sc]
        topClass = scDf["superclass"].value_counts().index[0] if not scDf.empty else ""
        net.add_node(
            f"SC_{sc}",
            label=label,
            title=f"[Scaffold] {sc}\n"
                  f"Compounds: {count}\n"
                  f"Superclass: {topClass}",
            color="#4a90d9",
            size=size,
            shape="dot",
            group="scaffold",
        )

    # 添加糖序列节点 (Add sugar nodes — medium, gold)
    sugarTotals = plotDf["sugar_sequence"].value_counts()
    for su in sugars:
        count = int(sugarTotals.get(su, 0))
        size = int(max(15, min(40, 10 + count * 2)))
        net.add_node(
            f"SU_{su}",
            label=su,
            title=f"[Sugar] {su}\nOccurrences: {count}",
            color="#f5a623",
            size=size,
            shape="dot",
            group="sugar",
        )

    # 添加物种节点 (Add organism nodes — small, green)
    orgTotals = plotDf["organism"].value_counts()
    for org in organisms:
        count = int(orgTotals.get(org, 0))
        size = int(max(8, min(25, 5 + count * 2)))
        # 只显示种加词 (Show species epithet for compact labels)
        parts = org.split()
        shortLabel = f"{parts[0][0]}. {parts[1]}" if len(parts) >= 2 else org
        net.add_node(
            f"ORG_{org}",
            label=shortLabel,
            title=f"[Organism] {org}\nCompounds: {count}",
            color="#2ecc71",
            size=size,
            shape="dot",
            group="organism",
        )

    # 添加边: Scaffold → Sugar (Add edges: Scaffold → Sugar)
    for (sc, su), count in scaffoldSugar.items():
        c = int(count)  # numpy int64 → Python int (pyvis JSON 兼容)
        width = float(max(1, min(6, c * 0.8)))
        net.add_edge(
            f"SC_{sc}", f"SU_{su}",
            value=c,
            width=width,
            color="rgba(74, 144, 217, 0.4)",
            title=f"{c} compounds",
        )

    # 添加边: Sugar → Organism (Add edges: Sugar → Organism)
    for (su, org), count in sugarOrganism.items():
        c = int(count)
        width = float(max(0.5, min(4, c * 0.5)))
        net.add_edge(
            f"SU_{su}", f"ORG_{org}",
            value=c,
            width=width,
            color="rgba(46, 204, 113, 0.3)",
            title=f"{c} compounds",
        )

    # 自定义 HTML 标题 (Custom HTML heading)
    net.heading = (
        f"<h2 style='color:#eee; text-align:center; font-family:Segoe UI,Arial;'>"
        f"{targetFamily} - Scaffold / Sugar / Organism Network</h2>"
        f"<p style='color:#aaa; text-align:center; font-size:13px;'>"
        f"<span style='color:#4a90d9'>&#9679;</span> Scaffold ({len(scaffolds)}) &nbsp;|&nbsp; "
        f"<span style='color:#f5a623'>&#9679;</span> Sugar Sequence ({len(sugars)}) &nbsp;|&nbsp; "
        f"<span style='color:#2ecc71'>&#9679;</span> Organism ({len(organisms)}) &nbsp;|&nbsp; "
        f"Total edges: {len(scaffoldSugar) + len(sugarOrganism)}</p>"
    )

    # 配置交互选项 (Configure interaction options)
    net.set_options("""
    {
        "interaction": {
            "hover": true,
            "navigationButtons": true,
            "tooltipDelay": 100,
            "zoomView": true
        },
        "physics": {
            "forceAtlas2Based": {
                "gravitationalConstant": -80,
                "centralGravity": 0.008,
                "springLength": 150,
                "springConstant": 0.04,
                "damping": 0.4
            },
            "minVelocity": 0.75,
            "solver": "forceAtlas2Based",
            "stabilization": { "iterations": 200 }
        }
    }
    """)

    # 手动 UTF-8 写入, 避免 Windows GBK 编码崩溃
    # Manual UTF-8 write to avoid Windows GBK crash
    htmlContent = net.generate_html()
    with open(outputPath, "w", encoding="utf-8") as f:
        f.write(htmlContent)
    print(f"  Network saved: {outputPath}")
    print(f"    Scaffolds: {len(scaffolds)}")
    print(f"    Sugar Sequences: {len(sugars)}")
    print(f"    Organisms: {len(organisms)}")
    print(f"    Edges: {len(scaffoldSugar) + len(sugarOrganism)}")

    return outputPath


# =====================================================================
# 4. 统计摘要 (Statistical Summary)
# =====================================================================

def printFamilySummary(familyDf: pd.DataFrame, targetFamily: str):
    """打印科级统计摘要 (Print family-level summary)."""
    print(f"\n{'='*60}")
    print(f"  {targetFamily} — Summary")
    print(f"{'='*60}")
    print(f"  Total compound-organism pairs: {len(familyDf)}")
    print(f"  Unique organisms: {familyDf['organism'].nunique()}")
    print(f"  Unique scaffolds: {familyDf['scaffold'].nunique()}")
    print(f"  Unique sugar sequences: {familyDf['sugar_sequence'].nunique()}")

    print(f"\n  [Top 5 Sugar Sequences]")
    for seq, count in familyDf["sugar_sequence"].value_counts().head(5).items():
        print(f"    {seq:40s} {count}")

    print(f"\n  [Top 5 Scaffolds]")
    for sc, count in familyDf["scaffold"].value_counts().head(5).items():
        label = sc[:50] + "..." if len(sc) > 50 else sc
        print(f"    {label:53s} {count}")

    print(f"\n  [Top 5 Organisms]")
    for org, count in familyDf["organism"].value_counts().head(5).items():
        print(f"    {org:40s} {count}")


# =====================================================================
# Main
# =====================================================================

def main():
    parser = argparse.ArgumentParser(description="GlycoNP Family Network Plot")
    parser.add_argument("--input", type=str, default=None, help="Input CSV")
    parser.add_argument("--family", type=str, default="Solanaceae", help="Target family")
    parser.add_argument("--min-scaffold", type=int, default=2,
                        help="Min compounds per scaffold to include")
    parser.add_argument("--all-families", action="store_true",
                        help="Generate for all major families")
    args = parser.parse_args()

    baseDir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    reportDir = os.path.join(baseDir, "reports")

    if args.input:
        inputPath = args.input
    else:
        inputPath = os.path.join(reportDir, "GlycoNP_Pipeline_Full.csv")
        if not os.path.exists(inputPath):
            inputPath = os.path.join(reportDir, "GlycoNP_Pipeline_1000.csv")

    print("=" * 60)
    print("  GlycoNP Family Network — Scaffold ↔ Sugar ↔ Organism")
    print("=" * 60)
    print(f"  Input: {inputPath}")

    df = pd.read_csv(inputPath, low_memory=False, dtype=str, encoding="utf-8-sig")
    print(f"  Loaded {len(df):,} rows")

    if args.all_families:
        families = [
            "Solanaceae", "Fabaceae", "Asteraceae", "Araliaceae",
            "Brassicaceae", "Rosaceae", "Lamiaceae", "Rutaceae",
        ]
    else:
        families = [args.family]

    for family in families:
        familyDf = filterByFamily(df, family)
        if familyDf.empty:
            print(f"\n  [SKIP] No compounds found for {family}")
            continue

        printFamilySummary(familyDf, family)

        outputPath = os.path.join(reportDir, f"Network_{family}.html")
        buildFamilyNetwork(familyDf, family, outputPath,
                           minScaffoldCount=args.min_scaffold)

    print(f"\n{'='*60}")
    print(f"  Done!")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
