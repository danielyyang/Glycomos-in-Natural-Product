"""
GlycoNP 修饰景观分析 — 同糖不同修
GlycoNP Modification Landscape — Same Sugar, Different Decorations

核心科学问题 (Core Scientific Question):
  "当糖链序列 (Sugar_Sequence) 相同时,
   大自然在不同骨架/物种/大类上给它加了什么不同的修饰?"

这揭示了:
  - 骨架驱动的修饰偏好 (Scaffold-driven modification preference)
  - 物种特异性的酶学指纹 (Species-specific enzymatic fingerprint)
  - 进化适应压力下的糖链"化妆"策略

使用方法 (Usage):
  python scripts/analyze_modification_landscape.py [--input PATH]
"""
import argparse
import os
import sys
from collections import Counter, defaultdict
from typing import Dict, List, Set, Tuple

import pandas as pd

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))


# 已知修饰基团列表 (Known modification types)
MOD_TYPES = ["O-Me", "O-Ac", "NAc", "Sulfate", "Phosphate", "COOH", "NH2"]


# =====================================================================
# 1. 修饰指纹提取 (Modification Fingerprint Extraction)
# =====================================================================

def parseModifications(modString: str) -> Dict[str, int]:
    """
    解析 Glycan_Modifications 字符串为修饰字典。
    Parse 'O-Me: 2, NAc: 1' into {'O-Me': 2, 'NAc': 1}.
    """
    result = {}
    if not modString or str(modString) in ("nan", "None", ""):
        return result

    for part in str(modString).split(","):
        part = part.strip()
        if ":" in part:
            key, val = part.rsplit(":", 1)
            key = key.strip()
            try:
                result[key] = int(val.strip())
            except ValueError:
                result[key] = 1
        elif part:
            # 没有计数的修饰 (Modification without count)
            result[part.strip()] = 1

    return result


def modFingerprint(modDict: Dict[str, int]) -> str:
    """
    将修饰字典转为标准化指纹字符串。
    Convert modification dict to normalized fingerprint string.
    {'O-Me': 2, 'NAc': 1} → 'NAc+O-Me×2'

    排序确保同一修饰组合生成相同指纹。
    """
    if not modDict:
        return "Unmodified"

    parts = []
    for mod in sorted(modDict.keys()):
        count = modDict[mod]
        if count > 1:
            parts.append(f"{mod}x{count}")
        else:
            parts.append(mod)

    return " + ".join(parts)


def cleanSuperclass(val: str) -> str:
    """清理 Superclass 名称。"""
    if not val or val == "nan":
        return "Unclassified"
    s = str(val).strip()
    if "(Tanimoto=" in s:
        s = s[:s.index("(Tanimoto=")].strip()
    if s.startswith("Glycolipid"):
        s = "Glycolipids"
    return s if s else "Unclassified"


def inferKingdom(organism: str) -> str:
    """从物种名推断生物界。"""
    if not organism or str(organism) in ("nan", ""):
        return "Unknown"
    org = str(organism).lower()
    bacteriaGenera = [
        "streptomyces", "bacillus", "pseudomonas", "escherichia",
        "staphylococcus", "mycobacterium", "salmonella", "lactobacillus",
    ]
    fungiGenera = [
        "aspergillus", "penicillium", "fusarium", "saccharomyces",
        "candida", "trichoderma", "agaricus", "ganoderma",
    ]
    marineKeywords = [
        "sponge", "coral", "halichondria", "aplysina", "ircinia",
        "haliclona", "dysidea", "axinella",
    ]
    for g in bacteriaGenera:
        if g in org:
            return "Bacteria"
    for g in fungiGenera:
        if g in org:
            return "Fungi"
    for kw in marineKeywords:
        if kw in org:
            return "Marine"
    return "Plantae"


# =====================================================================
# 2. 核心分析 (Core Analysis)
# =====================================================================

def analyzeSameSugarDiffMods(df: pd.DataFrame) -> List[str]:
    """
    核心: 同一 Sugar_Sequence 下, 修饰组合的多样性分析。
    Core: Modification diversity within the SAME Sugar_Sequence.

    对每条高频糖链, 统计:
      - 出现了几种不同的修饰 "指纹" (modification fingerprint)
      - 这些指纹分别在哪些大类/物种中出现
      - 哪些大类独占某种修饰?
    """
    # 只分析有意义的糖链 (Filter valid sequences)
    validDf = df[
        df["Sugar_Sequence"].notna()
        & (df["Sugar_Sequence"] != "")
        & (df["Sugar_Sequence"] != "nan")
    ].copy()

    validDf["_mods"] = validDf["Glycan_Modifications"].apply(
        lambda x: parseModifications(str(x))
    )
    validDf["_fp"] = validDf["_mods"].apply(modFingerprint)
    validDf["_class"] = validDf["Superclass"].apply(cleanSuperclass)

    # 获取 Top 15 高频糖链 (Top 15 sugar sequences by frequency)
    topSugars = validDf["Sugar_Sequence"].value_counts().head(15).index.tolist()

    lines = []
    lines.append("# Modification Landscape: Same Sugar, Different Decorations")
    lines.append("# 修饰景观: 同糖不同修\n")
    lines.append("> Core question: When the sugar sequence is identical,")
    lines.append("> how do modifications differ across scaffolds and species?\n")
    lines.append(f"> Total valid compounds: **{len(validDf):,}**\n")

    # ------------------------------------------------------------------
    # A. 每条糖链的修饰多样性 (Modification diversity per sugar)
    # ------------------------------------------------------------------
    lines.append("## A. Modification Diversity per Sugar Sequence")
    lines.append("## A. 每条糖链的修饰多样性\n")

    lines.append("| Sugar Sequence | Total | Unmodified | Modified | "
                 "Unique Mod Fingerprints | Diversity Index |")
    lines.append("|:--------------|------:|-----------:|---------:|"
                 "-----------------------:|----------------:|")

    diversityData = []

    for sugar in topSugars:
        subset = validDf[validDf["Sugar_Sequence"] == sugar]
        total = len(subset)
        unmod = len(subset[subset["_fp"] == "Unmodified"])
        mod = total - unmod
        uniqueFps = subset["_fp"].nunique()
        # Shannon 多样性指数 (Shannon diversity index)
        fpCounts = subset["_fp"].value_counts()
        probs = fpCounts / fpCounts.sum()
        shannon = -(probs * probs.apply(lambda p: p if p > 0 else 1e-10).apply(
            lambda p: __import__("math").log2(p)
        )).sum()

        diversityData.append({
            "sugar": sugar, "total": total, "unmod": unmod,
            "mod": mod, "uniqueFps": uniqueFps, "shannon": shannon,
        })

        lines.append(
            f"| `{sugar[:30]}` | {total:,} | {unmod:,} ({unmod/total*100:.0f}%) | "
            f"{mod:,} ({mod/total*100:.0f}%) | {uniqueFps} | "
            f"{shannon:.2f} |"
        )

    lines.append("")

    # ------------------------------------------------------------------
    # B. 同糖不同修 × 大类交叉 (Same sugar, diff mods × Superclass)
    # ------------------------------------------------------------------
    lines.append("## B. Modification Profiles by Superclass (Same Sugar)")
    lines.append("## B. 同一糖链在不同大类中的修饰画像\n")

    # 聚焦 Top 5 糖链 (Focus on top 5 sugars)
    for sugar in topSugars[:5]:
        subset = validDf[validDf["Sugar_Sequence"] == sugar]
        lines.append(f"### `{sugar}` (n={len(subset):,})\n")

        # 每个大类的修饰分布 (Mod distribution per superclass)
        topClasses = subset["_class"].value_counts().head(8).index.tolist()

        lines.append("| Superclass | n | Top Modification | "
                     "Modified % | Unique Fingerprints |")
        lines.append("|:-----------|--:|:----------------|----------:|"
                     "--------------------:|")

        for cls in topClasses:
            clsSubset = subset[subset["_class"] == cls]
            n = len(clsSubset)
            modSubset = clsSubset[clsSubset["_fp"] != "Unmodified"]
            modPct = len(modSubset) / n * 100 if n > 0 else 0

            # 该大类最常见的修饰 (Most common mod in this class)
            allMods = Counter()
            for mods in clsSubset["_mods"]:
                for mod, count in mods.items():
                    allMods[mod] += count
            topMod = allMods.most_common(1)
            topModStr = f"{topMod[0][0]} ({topMod[0][1]})" if topMod else "None"

            uniqueFps = clsSubset["_fp"].nunique()
            lines.append(
                f"| {cls[:25]} | {n} | {topModStr} | "
                f"{modPct:.0f}% | {uniqueFps} |"
            )

        lines.append("")

    # ------------------------------------------------------------------
    # C. 同糖不同修 × 物种界 (Same sugar, diff mods × Kingdom)
    # ------------------------------------------------------------------
    lines.append("## C. Modification Profiles by Kingdom (Same Sugar)")
    lines.append("## C. 同一糖链在不同生物界中的修饰差异\n")

    # 展开物种行 (Explode organisms)
    orgRows = []
    for _, row in validDf.iterrows():
        orgs = str(row.get("organisms", "")).split("|")
        for org in orgs:
            org = org.strip()
            if not org or org == "nan":
                continue
            orgRows.append({
                "Sugar_Sequence": row["Sugar_Sequence"],
                "_fp": row["_fp"],
                "_mods": row["_mods"],
                "kingdom": inferKingdom(org),
            })

    if orgRows:
        orgDf = pd.DataFrame(orgRows)
        kingdoms = ["Plantae", "Bacteria", "Fungi", "Marine"]

        for sugar in topSugars[:5]:
            subset = orgDf[orgDf["Sugar_Sequence"] == sugar]
            if len(subset) < 10:
                continue

            lines.append(f"### `{sugar}` — Kingdom Comparison\n")
            lines.append("| Kingdom | n | Modified % | Top Mod | "
                         "Exclusive Mods |")
            lines.append("|:--------|--:|----------:|:--------|"
                         ":--------------|")

            # 收集每个界的修饰集合 (Collect mod sets per kingdom)
            kingdomMods: Dict[str, Set[str]] = {}

            for kingdom in kingdoms:
                kSubset = subset[subset["kingdom"] == kingdom]
                if kSubset.empty:
                    continue
                n = len(kSubset)
                modPct = len(kSubset[kSubset["_fp"] != "Unmodified"]) / n * 100

                allMods = Counter()
                modSet = set()
                for mods in kSubset["_mods"]:
                    for mod, cnt in mods.items():
                        allMods[mod] += cnt
                        modSet.add(mod)

                kingdomMods[kingdom] = modSet
                topMod = allMods.most_common(1)
                topModStr = topMod[0][0] if topMod else "None"

                lines.append(
                    f"| {kingdom} | {n} | {modPct:.0f}% | {topModStr} | - |"
                )

            # 回填独占修饰 (Fill in exclusive mods)
            # 那些只在某个界出现, 其他界没有的修饰
            allKingdomMods = set()
            for mSet in kingdomMods.values():
                allKingdomMods |= mSet

            exclusives = {}
            for kingdom, mSet in kingdomMods.items():
                otherMods = set()
                for otherK, otherSet in kingdomMods.items():
                    if otherK != kingdom:
                        otherMods |= otherSet
                excl = mSet - otherMods
                if excl:
                    exclusives[kingdom] = ", ".join(sorted(excl))

            if exclusives:
                lines.append(f"\n**Exclusive modifications** (only in one kingdom):")
                for k, mods in exclusives.items():
                    lines.append(f"  - **{k}**: {mods}")

            lines.append("")

    # ------------------------------------------------------------------
    # D. 修饰组合热点 (Modification Combination Hotspots)
    # ------------------------------------------------------------------
    lines.append("## D. Modification Combination Hotspots")
    lines.append("## D. 修饰组合热点\n")
    lines.append("Which modification combinations most frequently decorate "
                 "the same sugar?\n")

    # 过滤有修饰的行 (Filter modified rows)
    modOnly = validDf[validDf["_fp"] != "Unmodified"].copy()
    if not modOnly.empty:
        fpCounts = modOnly["_fp"].value_counts().head(15)
        lines.append("| Rank | Modification Fingerprint | Count | "
                     "Top Sugar | Top Superclass |")
        lines.append("|:----:|:------------------------|------:|"
                     ":----------|:---------------|")

        for i, (fp, count) in enumerate(fpCounts.items(), 1):
            fpSubset = modOnly[modOnly["_fp"] == fp]
            topSugar = fpSubset["Sugar_Sequence"].value_counts().index[0]
            topClass = fpSubset["_class"].value_counts().index[0]
            lines.append(
                f"| {i} | `{fp}` | {count:,} | "
                f"`{topSugar[:25]}` | {topClass[:20]} |"
            )
        lines.append("")

    # ------------------------------------------------------------------
    # E. 进化密码 (Evolutionary Code)
    # ------------------------------------------------------------------
    lines.append("## E. Evolutionary Insights / 进化密码\n")
    lines.append("Sugar sequences where modifications differ MOST between "
                 "superclasses:\n")

    # 找修饰指纹在不同大类间差异最大的糖链
    # (Find sugars where mod fingerprints vary most across superclasses)
    candidateRows = []
    for sugar in topSugars[:10]:
        subset = validDf[validDf["Sugar_Sequence"] == sugar]
        if len(subset) < 20:
            continue

        # 每个大类的修饰率 (Mod rate per class)
        classGroups = subset.groupby("_class")
        classModRates = {}
        for cls, grp in classGroups:
            if len(grp) < 5:
                continue
            modRate = len(grp[grp["_fp"] != "Unmodified"]) / len(grp)
            classModRates[cls] = modRate

        if len(classModRates) >= 2:
            rates = list(classModRates.values())
            maxDiff = max(rates) - min(rates)
            highCls = max(classModRates, key=classModRates.get)
            lowCls = min(classModRates, key=classModRates.get)
            candidateRows.append({
                "sugar": sugar,
                "maxDiff": maxDiff,
                "highClass": highCls,
                "highRate": classModRates[highCls],
                "lowClass": lowCls,
                "lowRate": classModRates[lowCls],
            })

    if candidateRows:
        candidateRows.sort(key=lambda x: -x["maxDiff"])

        lines.append("| Sugar Sequence | Most Modified Class | % | "
                     "Least Modified Class | % | Delta |")
        lines.append("|:--------------|:-------------------|--:|"
                     ":--------------------|--:|------:|")

        for row in candidateRows[:8]:
            lines.append(
                f"| `{row['sugar'][:25]}` | {row['highClass'][:20]} | "
                f"{row['highRate']*100:.0f}% | {row['lowClass'][:20]} | "
                f"{row['lowRate']*100:.0f}% | "
                f"{row['maxDiff']*100:.0f}pp |"
            )
        lines.append("")

    lines.append("---")
    lines.append("> Generated by `analyze_modification_landscape.py` | "
                 "GlycoNP Pipeline")

    return lines


# =====================================================================
# Main
# =====================================================================

def main():
    parser = argparse.ArgumentParser(
        description="GlycoNP Modification Landscape Analysis")
    parser.add_argument("--input", type=str, default=None, help="Input CSV")
    args = parser.parse_args()

    baseDir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    reportDir = os.path.join(baseDir, "reports")

    if args.input:
        inputPath = args.input
    else:
        inputPath = os.path.join(reportDir, "GlycoNP_Pipeline_Full.csv")
        if not os.path.exists(inputPath):
            inputPath = os.path.join(reportDir, "GlycoNP_Pipeline_1000.csv")

    print("=" * 70)
    print("  GlycoNP Modification Landscape Analysis")
    print("  同糖不同修 — Same Sugar, Different Decorations")
    print("=" * 70)
    print(f"  Input: {inputPath}")

    df = pd.read_csv(inputPath, low_memory=False, dtype=str, encoding="utf-8-sig")
    print(f"  Loaded {len(df):,} rows")

    lines = analyzeSameSugarDiffMods(df)

    reportPath = os.path.join(reportDir, "Modification_Landscape.md")
    with open(reportPath, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))

    print(f"\n  Report saved: {reportPath}")
    print(f"  Total lines: {len(lines)}")
    print(f"\n{'='*70}")
    print(f"  Done!")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
