"""
Phase 4 Family Gap-Fill — 基于化合物名称的科级填补
Phase 4 Family Gap-Fill — Name-Based Family Assignment

针对 30.6% 缺失 Family 的化合物, 通过名称关键词匹配推断所属科。
For compounds missing Family data, infer family from compound name keywords.

策略 (Strategy):
  1. 内建名称→科映射字典 (100+ 常见天然产物根词)
  2. 按 name 列进行关键词匹配
  3. 可选: pubchempy 二级查询

使用方法 (Usage):
  python scripts/fill_family_by_name.py [--input PATH]
"""
import argparse
import os
import re
import sys
from typing import Dict, Optional

import pandas as pd

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))


# =====================================================================
# 名称→科 映射字典 (Name → Family dictionary)
# 基于天然产物化学领域的经典命名规则
# =====================================================================

NAME_TO_FAMILY: Dict[str, str] = {
    # 人参/五加科 (Araliaceae — Ginseng family)
    "ginsenoside": "Araliaceae",
    "panax": "Araliaceae",
    "araloside": "Araliaceae",
    "eleutherosi": "Araliaceae",
    "chikusetsusaponin": "Araliaceae",

    # 茄科 (Solanaceae — Nightshade)
    "solanine": "Solanaceae",
    "tomatine": "Solanaceae",
    "solasonine": "Solanaceae",
    "solasodine": "Solanaceae",
    "chaconine": "Solanaceae",
    "capsaicin": "Solanaceae",
    "withanolide": "Solanaceae",
    "datura": "Solanaceae",
    "nicotine": "Solanaceae",

    # 豆科 (Fabaceae — Legume)
    "soyasaponin": "Fabaceae",
    "glycyrrhiz": "Fabaceae",
    "liquirit": "Fabaceae",
    "astragaloside": "Fabaceae",
    "isoflavone": "Fabaceae",
    "geniste": "Fabaceae",
    "daidzein": "Fabaceae",
    "formononetin": "Fabaceae",
    "trifolirhizin": "Fabaceae",
    "sophoraflavanone": "Fabaceae",

    # 十字花科 (Brassicaceae — Crucifer)
    "glucosinolate": "Brassicaceae",
    "sinigrin": "Brassicaceae",
    "glucoraphanin": "Brassicaceae",
    "sulforaphane": "Brassicaceae",

    # 菊科 (Asteraceae — Composite)
    "artemisinin": "Asteraceae",
    "parthenolide": "Asteraceae",
    "helianth": "Asteraceae",
    "inulin": "Asteraceae",
    "taraxast": "Asteraceae",
    "lactuci": "Asteraceae",
    "sesquiterpene lactone": "Asteraceae",
    "cynarin": "Asteraceae",

    # 芸香科 (Rutaceae — Citrus)
    "hesperidin": "Rutaceae",
    "naringin": "Rutaceae",
    "nobiletin": "Rutaceae",
    "limonin": "Rutaceae",
    "citrus": "Rutaceae",
    "tangeretin": "Rutaceae",
    "aurantin": "Rutaceae",

    # 蔷薇科 (Rosaceae — Rose)
    "amygdalin": "Rosaceae",
    "prunin": "Rosaceae",
    "rosarin": "Rosaceae",
    "prunasin": "Rosaceae",
    "tormentic": "Rosaceae",

    # 唇形科 (Lamiaceae — Mint)
    "rosmarinic": "Lamiaceae",
    "salvianolic": "Lamiaceae",
    "tanshinone": "Lamiaceae",
    "baicalin": "Lamiaceae",
    "scutellarin": "Lamiaceae",
    "acteoside": "Lamiaceae",

    # 百合科/石蒜科 (Amaryllidaceae)
    "galantamine": "Amaryllidaceae",
    "lycorine": "Amaryllidaceae",

    # 伞形科 (Apiaceae)
    "saponin": "Apiaceae",
    "angelicin": "Apiaceae",
    "imperatorin": "Apiaceae",

    # 毛茛科 (Ranunculaceae)
    "berberine": "Ranunculaceae",
    "coptisine": "Ranunculaceae",
    "aconitine": "Ranunculaceae",

    # 木犀科 (Oleaceae)
    "oleuropein": "Oleaceae",
    "ligustroside": "Oleaceae",
    "oleuroside": "Oleaceae",

    # 禾本科 (Poaceae — Grass)
    "oryzanol": "Poaceae",
    "tricin": "Poaceae",

    # 葡萄科 (Vitaceae)
    "resveratrol": "Vitaceae",
    "viniferin": "Vitaceae",

    # 大戟科 (Euphorbiaceae)
    "phorbol": "Euphorbiaceae",
    "euphorbia": "Euphorbiaceae",
    "ricin": "Euphorbiaceae",

    # 番荔枝科 (Annonaceae)
    "annonacin": "Annonaceae",
    "squamocin": "Annonaceae",

    # 龙胆科 (Gentianaceae)
    "gentiopicroside": "Gentianaceae",
    "swertiamarin": "Gentianaceae",

    # 鼠李科 (Rhamnaceae)
    "jujuboside": "Rhamnaceae",
    "ziziphin": "Rhamnaceae",

    # 茜草科 (Rubiaceae)
    "aucubin": "Rubiaceae",
    "quinine": "Rubiaceae",
    "catharanthine": "Rubiaceae",

    # 细菌/微生物 (Bacteria — not a plant family, but useful marker)
    "streptomycin": "Streptomycetaceae",
    "kanamycin": "Streptomycetaceae",
    "erythromycin": "Streptomycetaceae",
    "vancomycin": "Streptomycetaceae",
    "daunorubicin": "Streptomycetaceae",
}


def matchFamilyByName(name: str) -> Optional[str]:
    """
    通过化合物名称关键词匹配推断科名。
    Infer family from compound name via keyword matching.
    """
    if not name or str(name) in ("nan", ""):
        return None

    nameLower = str(name).lower().strip()

    for keyword, family in NAME_TO_FAMILY.items():
        if keyword in nameLower:
            return family

    return None


def main():
    parser = argparse.ArgumentParser(
        description="Phase 4 Family Gap-Fill by Name")
    parser.add_argument("--input", type=str, default=None, help="Input CSV")
    parser.add_argument("--dry-run", action="store_true",
                        help="Preview only, don't write")
    args = parser.parse_args()

    baseDir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    reportDir = os.path.join(baseDir, "reports")

    if args.input:
        inputPath = args.input
    else:
        inputPath = os.path.join(reportDir, "GlycoNP_Pipeline_Full.csv")

    print("=" * 60)
    print("  Phase 4 Family Gap-Fill by Name")
    print("  基于名称的科级填补")
    print("=" * 60)
    print(f"  Input: {inputPath}")

    df = pd.read_csv(inputPath, low_memory=False, dtype=str,
                     encoding="utf-8-sig")
    total = len(df)
    familyCol = "organism_taxonomy_05family"

    # 如果列不存在则创建 (Create column if missing)
    if familyCol not in df.columns:
        df[familyCol] = ""

    # 统计填补前 (Before stats)
    missingBefore = df[familyCol].isna() | (df[familyCol] == "") | (df[familyCol] == "nan")
    missingCountBefore = missingBefore.sum()
    coverageBefore = (total - missingCountBefore) / total * 100

    print(f"  Total: {total:,}")
    print(f"  Missing family BEFORE: {missingCountBefore:,} ({100-coverageBefore:.1f}%)")
    print(f"  Coverage BEFORE: {coverageBefore:.1f}%")

    # 执行名称匹配 (Run name matching)
    filled = 0
    familyStats = {}
    for idx in df[missingBefore].index:
        name = str(df.at[idx, "name"])
        family = matchFamilyByName(name)
        if family:
            df.at[idx, familyCol] = family
            filled += 1
            familyStats[family] = familyStats.get(family, 0) + 1

    # 也尝试 iupac_name (Try IUPAC name too)
    missingAfterName = df[familyCol].isna() | (df[familyCol] == "") | (df[familyCol] == "nan")
    for idx in df[missingAfterName].index:
        iupac = str(df.at[idx, "iupac_name"]) if "iupac_name" in df.columns else ""
        family = matchFamilyByName(iupac)
        if family:
            df.at[idx, familyCol] = family
            filled += 1
            familyStats[family] = familyStats.get(family, 0) + 1

    # 统计填补后 (After stats)
    missingAfter = df[familyCol].isna() | (df[familyCol] == "") | (df[familyCol] == "nan")
    missingCountAfter = missingAfter.sum()
    coverageAfter = (total - missingCountAfter) / total * 100

    print(f"\n  Filled by name: {filled:,}")
    print(f"  Missing family AFTER: {missingCountAfter:,} ({100-coverageAfter:.1f}%)")
    print(f"  Coverage AFTER: {coverageAfter:.1f}%")
    print(f"  Improvement: +{coverageAfter - coverageBefore:.1f}pp")

    if familyStats:
        print(f"\n  Top fills:")
        for fam, cnt in sorted(familyStats.items(), key=lambda x: -x[1])[:10]:
            print(f"    {fam:30s} {cnt}")

    if not args.dry_run:
        # 输出到新文件, 避免锁定冲突 (Output to new file to avoid lock)
        outPath = inputPath.replace(".csv", "_GapFilled.csv")
        df.to_csv(outPath, index=False, encoding="utf-8-sig")
        print(f"\n  Saved: {outPath}")
    else:
        print(f"\n  [DRY RUN] No changes written")

    print(f"\n{'='*60}")


if __name__ == "__main__":
    main()
