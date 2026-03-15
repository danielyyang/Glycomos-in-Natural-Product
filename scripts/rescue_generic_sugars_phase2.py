"""
GlycoNP Phase 2 深度泛指糖修复 (Deep Generic Sugar Rescue)

针对 Phase 1 未修复的 ~31K 泛指糖 tokens, 采用三种降级策略:
  E. IUPAC 立体化学深度解析 (tetrahydropyran + chiral centers)
  F. 高碳糖生物合成法则 (Non→Neu5Ac, Oct→KDO)
  G. InChIKey First Block 二维骨架跨库关联

使用方法 (Usage):
  python scripts/rescue_generic_sugars_phase2.py
"""
import os
import re
import sys
import time
from collections import Counter

import pandas as pd
from rdkit import Chem
from rdkit.Chem.inchi import MolToInchi, InchiToInchiKey

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

# =====================================================================
# 策略 E: IUPAC 立体化学手性序列解析
# =====================================================================
# 经典吡喃单糖的 IUPAC 手性指纹 (Canonical chiral fingerprints)
# 格式: 糖环碳 C2-C5 的 R/S 序列 (在 tetrahydropyran 上下文中)
#
# 化学原理:
#   D-Glucose pyranose:  C2=R, C3=S, C4=S, C5=R  (最常见的 IUPAC 编码)
#   D-Galactose pyranose: C2=R, C3=S, C4=R, C5=R  (C4 差向)
#   D-Mannose pyranose:  C2=S, C3=S, C4=S, C5=R  (C2 差向)
#   D-Xylose pyranose:   C2=R, C3=S, C4=S        (五碳, 无 C6)
#   L-Rhamnose pyranose: C2=S, C3=R, C4=R, C5=S  (6-脱氧-L-甘露糖)
#   L-Arabinose pyranose: C2=S, C3=R, C4=R       (五碳)
#   L-Fucose pyranose:   C2=S, C3=R, C4=S, C5=S  (6-脱氧-L-半乳糖)
#
# 注意: IUPAC 命名中的编号与 Fischer 投影可能不完全一致,
#       但 "trihydroxyoxan" 模式下相对构型是可靠的。

# 正则: 提取 tetrahydropyran 附近的手性序列
# 支持两种格式: (2R,3S,4S,5R) 和 (2~{R},3~{S},4~{S},5~{R})
CHIRAL_RE = re.compile(
    r'\('
    r'(\d)~?\{?([RS])\}?'
    r'(?:,(\d)~?\{?([RS])\}?)*'
    r'\)'
)

# 更健壮的提取: 找到所有 nR/nS 对
STEREOCENTER_RE = re.compile(r'(\d)~?\{?([RS])\}?')

# 吡喃糖的标志词 (tetrahydropyran variants in IUPAC)
PYRANOSE_INDICATORS = [
    'tetrahydropyran', 'oxan-2-', 'oxane',
    'trihydroxyoxan', 'trihydroxy-6-',
    'tetrahydro-2h-pyran',
]
FURANOSE_INDICATORS = [
    'tetrahydrofuran', 'oxolan', 'oxolane',
]

# 手性指纹字典 (chiral fingerprint → sugar name)
# 键: 冻结集合 {(position, chirality)} 的 sorted tuple
# 从 trihydroxyoxan (三羟基吡喃) 上下文中提取
#
# 六碳吡喃糖: C2-C5 四个手性中心
PYRANOSE_CHIRAL_MAP = {
    # D-Glc: 经典的 (3,4,5-trihydroxy-6-(hydroxymethyl)tetrahydropyran-2-yl)
    # C位模式: 2R,3R,4S,5S,6R 或 2R,4S,5R (简化)
    ("R", "R", "S", "S"): "D-Glc",    # C2R,C3R,C4S,C5S
    ("R", "S", "S", "R"): "D-Glc",    # 另一种编号
    ("S", "S", "R", "R"): "D-Glc",    # 镜像编号变体

    # D-Gal: C4 差向异构
    ("R", "R", "R", "S"): "D-Gal",    # C4 从 S→R
    ("R", "S", "R", "R"): "D-Gal",
    ("S", "R", "R", "R"): "D-Gal",

    # D-Man: C2 差向异构
    ("S", "R", "S", "S"): "D-Man",
    ("S", "S", "S", "R"): "D-Man",

    # L-Rha: 6-deoxy-L-mannose 特征
    ("S", "R", "R", "S"): "L-Rha",
    ("R", "R", "S", "S"): "L-Rha",     # 在 methyl 上下文中

    # D-Xyl: 五碳 (3个手性中心)
    ("R", "S", "S"): "D-Xyl",
    ("S", "R", "R"): "D-Xyl",

    # L-Ara: 五碳
    ("S", "R", "R"): "L-Ara",
    ("R", "S", "S"): "L-Ara",   # 上下文区分
}


def strategyE_iupacChiralParsing(
    iupacName: str,
    genericToken: str,
) -> str:
    """
    策略 E: 从 IUPAC 名中解析手性序列, 映射到精确单糖。
    Strategy E: Parse chiral sequences from IUPAC name for sugar ID.
    """
    if not iupacName or str(iupacName) == "nan":
        return None

    text = str(iupacName).lower()

    # 检查是否含有吡喃/呋喃描述
    hasPyranose = any(ind in text for ind in PYRANOSE_INDICATORS)
    hasFuranose = any(ind in text for ind in FURANOSE_INDICATORS)

    if not hasPyranose and not hasFuranose:
        return None

    # 快速规则: 最常见的 IUPAC 模式直接匹配
    # 模式 1: "3,4,5-trihydroxy-6-(hydroxymethyl)tetrahydropyran-2-yl"
    # 这是 D-Glc 的最典型 IUPAC 碎片
    if "trihydroxy-6-(hydroxymethyl)" in text and hasPyranose:
        # 提取附近的手性标记
        # 在 D-Glc 上下文中, 常见: (2R,3R,4S,5S,6R) 或 (2R,4S,5R)
        # 在 D-Gal 上下文中, 常见: (2R,3R,4R,5S,6R)
        # 检查 C4 的手性 → 区分 Glc (S) vs Gal (R)
        if re.search(r'4~?\{?S\}?.*trihydroxy|trihydroxy.*4~?\{?S\}?', text):
            return "D-Glc"
        elif re.search(r'4~?\{?R\}?.*trihydroxy|trihydroxy.*4~?\{?R\}?', text):
            return "D-Gal"
        else:
            # 无手性标记但有 trihydroxy-hydroxymethyl-pyran 模式
            # 统计上 D-Glc 占绝大多数
            return "D-Glc_iupac_parsed"

    # 模式 2: "3,4,5-trihydroxy-6-methyltetrahydropyran" (甲基=脱氧糖)
    if re.search(r'trihydroxy.*6-methyl.*(?:tetrahydropyran|oxan)', text):
        # 6-methyl = 6-deoxy 脱氧糖: L-Rha 或 L-Fuc
        if re.search(r'4~?\{?S\}?', text):
            return "L-Fuc"
        elif re.search(r'4~?\{?R\}?', text):
            return "L-Rha"
        else:
            return "L-Rha_iupac_parsed"  # 脱氧糖中 L-Rha 最常见

    # 模式 3: "3,4-dihydroxy-5-(hydroxymethyl)tetrahydrofuran" (呋喃糖)
    if hasFuranose and "dihydroxy" in text and "hydroxymethyl" in text:
        return "D-Ribf_iupac_parsed"  # 呋喃核糖最常见

    # 模式 4: "3,4,5-trihydroxytetrahydropyran" (无 hydroxymethyl = 五碳糖)
    if re.search(r'trihydroxy.*(?:tetrahydropyran|oxan)', text) and \
       "hydroxymethyl" not in text:
        if genericToken == "Pen":
            return "D-Xyl_iupac_parsed"  # 五碳吡喃糖最常见
        return "Pen_iupac_parsed"

    # 模式 5: "dihydroxy-6-(hydroxymethyl)tetrahydropyran" (二羟基六碳=脱氧糖)
    if re.search(r'dihydroxy.*6-\(hydroxymethyl\).*(?:tetrahydropyran|oxan)', text):
        return "dHex_iupac_parsed"

    return None


# =====================================================================
# 策略 F: 高碳糖生物合成法则
# =====================================================================
# Non (壬糖, 9C) + 特定大类 → Neu5Ac (唾液酸)
# Oct (辛糖, 8C) + 特定大类 → KDO (3-脱氧-D-甘露-辛-2-酮糖酸)

SIALIC_ACID_CLASSES = {
    "Glycolipid", "Sphingolipid", "Sphingoid_base",
    "Glycopeptide", "Amino acid glycosides",
}

def strategyF_biologicalInference(
    genericToken: str,
    superclass: str,
) -> str:
    """
    策略 F: 基于生物合成逻辑的强硬编码推断。
    Strategy F: Biological inference for high-carbon sugars.
    """
    scClean = str(superclass).strip() if superclass else ""
    if "(Tanimoto=" in scClean:
        scClean = scClean[:scClean.index("(Tanimoto=")].strip()

    # 规则 1: Non + 鞘脂/糖脂/糖肽 → Neu5Ac
    if genericToken == "Non":
        if any(key in scClean for key in SIALIC_ACID_CLASSES):
            return "Neu5Ac_biological_inferred"

    # 规则 2: Oct + 糖脂/细菌源 → KDO
    if genericToken == "Oct":
        if any(key in scClean for key in {"Glycolipid", "Saccharides"}):
            return "KDO_biological_inferred"

    # 规则 3: Non + 任何糖脂子类 (Aliphatic C*) → Neu5Ac
    if genericToken == "Non" and "Glycolipid" in scClean:
        return "Neu5Ac_biological_inferred"

    # 规则 4: Hept + Macrolides → 大环内酯庚糖
    if genericToken == "Hept" and "Macrolides" in scClean:
        return "L-Mycarose_biological_inferred"

    return None


# =====================================================================
# 策略 G: InChIKey First Block 二维骨架关联
# =====================================================================
def buildInchiKeyPrior(df: pd.DataFrame) -> dict:
    """
    构建 InChIKey First Block → 精确糖 的映射表。
    从已精确匹配的样本中, 按 InChIKey 14-char block 分组统计。
    """
    print("  [Strategy G] Building InChIKey first-block prior...")

    genericPat = r'\bHex\b|\bPen\b|\bdHex\b|\bNon\b|\bOct\b|\bHept\b'
    preciseMask = (
        df["Sugar_Sequence"].notna()
        & (df["Sugar_Sequence"] != "")
        & (~df["Sugar_Sequence"].str.contains(genericPat, na=False, regex=True))
        & (~df["Sugar_Sequence"].str.contains("_predicted|_iupac|_inferred",
                                               na=False, regex=True))
        & df["standard_inchi_key"].notna()
    )
    preciseDf = df.loc[preciseMask].copy()
    print(f"    Precise + InChIKey samples: {len(preciseDf):,}")

    def extractFirst(seq):
        tokens = re.findall(
            r'Neu5Ac|Neu5Gc|KDO|[DL]-[A-Z][a-z]+[A-Z]?[a-z]*',
            str(seq))
        return tokens[0] if tokens else None

    preciseDf["_sugar"] = preciseDf["Sugar_Sequence"].apply(extractFirst)
    preciseDf = preciseDf[preciseDf["_sugar"].notna()]

    # 截取 InChIKey First Block (前14字符 = 二维连接性)
    preciseDf["_fb"] = preciseDf["standard_inchi_key"].apply(
        lambda x: str(x)[:14] if pd.notna(x) and len(str(x)) >= 14 else None)
    preciseDf = preciseDf[preciseDf["_fb"].notna()]

    priorTable = {}
    for fb, group in preciseDf.groupby("_fb"):
        dist = group["_sugar"].value_counts().to_dict()
        priorTable[fb] = dist

    print(f"    First-block priors: {len(priorTable):,}")
    return priorTable


def strategyG_inchikeyMatch(
    smiles: str,
    inchiKey: str,
    priorTable: dict,
    minConfidence: float = 0.80,
) -> str:
    """
    策略 G: 通过 InChIKey First Block 跨库关联推断精确糖。
    Strategy G: InChIKey first-block 2D topology match.
    """
    fb = None

    # 优先使用已有的 InChIKey
    if inchiKey and str(inchiKey) != "nan" and len(str(inchiKey)) >= 14:
        fb = str(inchiKey)[:14]

    # 回退: 从 SMILES 计算
    if fb is None and smiles and str(smiles) != "nan":
        try:
            mol = Chem.MolFromSmiles(str(smiles))
            if mol:
                inchi = MolToInchi(mol)
                if inchi:
                    ik = InchiToInchiKey(inchi)
                    if ik and len(ik) >= 14:
                        fb = ik[:14]
        except Exception:
            pass

    if fb is None:
        return None

    dist = priorTable.get(fb, {})
    total = sum(dist.values())
    if total < 3:
        return None

    topSugar = max(dist, key=dist.get)
    confidence = dist[topSugar] / total
    if confidence >= minConfidence:
        return f"{topSugar}_2d_inferred"

    return None


# =====================================================================
# 主管线 (Main Pipeline)
# =====================================================================
GENERIC_PATTERN = re.compile(r'\b(Hex|Pen|dHex|HexA|Non|Oct|Hept)\b')

# 排除已修复的 token (带有 _predicted/_iupac/_inferred 后缀的不再处理)
ALREADY_RESCUED = re.compile(
    r'_predicted|_iupac_parsed|_biological_inferred|_2d_inferred')


def rescueTokenPhase2(
    token: str,
    iupacName: str,
    superclass: str,
    smiles: str,
    inchiKey: str,
    inchiKeyPrior: dict,
    logCounter: dict,
) -> str:
    """Phase 2 修复单个 token: E → F → G 级联"""
    if token not in ("Hex", "Pen", "dHex", "HexA", "Non", "Oct", "Hept"):
        return token

    # 策略 E: IUPAC 手性解析
    resultE = strategyE_iupacChiralParsing(iupacName, token)
    if resultE:
        logCounter["E"] += 1
        return resultE

    # 策略 F: 生物合成推断
    resultF = strategyF_biologicalInference(token, superclass)
    if resultF:
        logCounter["F"] += 1
        return resultF

    # 策略 G: InChIKey 2D 关联
    resultG = strategyG_inchikeyMatch(smiles, inchiKey, inchiKeyPrior)
    if resultG:
        logCounter["G"] += 1
        return resultG

    logCounter["MISS"] += 1
    return token


def rescueSequencePhase2(
    seq: str,
    iupacName: str,
    superclass: str,
    smiles: str,
    inchiKey: str,
    inchiKeyPrior: dict,
    logCounter: dict,
) -> str:
    """修复整条序列中的泛化 token"""
    def replaceMatch(m):
        return rescueTokenPhase2(
            m.group(0), iupacName, superclass, smiles, inchiKey,
            inchiKeyPrior, logCounter)
    return GENERIC_PATTERN.sub(replaceMatch, seq)


def main():
    baseDir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    reportDir = os.path.join(baseDir, "reports")
    inputPath = os.path.join(reportDir, "GlycoNP_Pipeline_Full_Cleaned.csv")

    print("=" * 70)
    print("  GlycoNP Phase 2 — Deep Generic Sugar Rescue")
    print("  深度泛指糖修复 (策略 E/F/G)")
    print("=" * 70)

    t0 = time.time()
    df = pd.read_csv(inputPath, low_memory=False, encoding="utf-8-sig")
    total = len(df)
    print(f"  Loaded: {total:,} rows ({time.time()-t0:.1f}s)")

    # 找到仍有泛化标签的行
    genericPat = r'\bHex\b|\bPen\b|\bdHex\b|\bHexA\b|\bNon\b|\bOct\b|\bHept\b'
    targetMask = df["Sugar_Sequence"].str.contains(
        genericPat, na=False, regex=True)
    targetCount = targetMask.sum()

    # 旧统计
    oldTokens = []
    for seq in df.loc[targetMask, "Sugar_Sequence"].dropna():
        oldTokens.extend(re.findall(genericPat, str(seq)))
    oldCounter = Counter(oldTokens)
    totalOld = sum(oldCounter.values())

    print(f"\n  [AUDIT] Remaining generic tokens: {totalOld:,} in {targetCount:,} rows")
    for k, v in oldCounter.most_common():
        print(f"    {k:10s} {v:>6,}")

    # 构建策略 G 先验
    inchiKeyPrior = buildInchiKeyPrior(df)

    # 执行修复
    print(f"\n  [EXECUTE] Phase 2 Rescue: E → F → G cascade...")
    logCounter = {"E": 0, "F": 0, "G": 0, "MISS": 0}
    rescued = 0

    for idx in df.index[targetMask]:
        oldSeq = str(df.at[idx, "Sugar_Sequence"])
        newSeq = rescueSequencePhase2(
            oldSeq,
            str(df.at[idx, "iupac_name"]) if pd.notna(df.at[idx, "iupac_name"]) else None,
            str(df.at[idx, "Superclass"]) if pd.notna(df.at[idx, "Superclass"]) else None,
            str(df.at[idx, "canonical_smiles"]) if pd.notna(df.at[idx, "canonical_smiles"]) else None,
            str(df.at[idx, "standard_inchi_key"]) if pd.notna(df.at[idx, "standard_inchi_key"]) else None,
            inchiKeyPrior, logCounter,
        )
        if newSeq != oldSeq:
            df.at[idx, "Sugar_Sequence"] = newSeq
            rescued += 1

    # 结果汇总
    totalRescued = logCounter["E"] + logCounter["F"] + logCounter["G"]
    rescueRate = totalRescued / totalOld * 100 if totalOld > 0 else 0

    print(f"\n  [RESULT] Phase 2 Rescue Summary:")
    print(f"    Rows modified: {rescued:,} / {targetCount:,}")
    print(f"    Strategy E (IUPAC Chiral):      {logCounter['E']:>6,} tokens")
    print(f"    Strategy F (Bio-Inference):      {logCounter['F']:>6,} tokens")
    print(f"    Strategy G (InChIKey 2D):        {logCounter['G']:>6,} tokens")
    print(f"    Unrescued (MISS):                {logCounter['MISS']:>6,} tokens")
    print(f"    Phase 2 rescue rate: {rescueRate:.1f}%")

    # 新统计
    newTokens = []
    for seq in df["Sugar_Sequence"].dropna():
        newTokens.extend(re.findall(genericPat, str(seq)))
    newCounter = Counter(newTokens)

    print(f"\n  [COMPARISON] Before → After:")
    print(f"  {'Token':<12s} {'Before':>8s} {'After':>8s} {'Delta':>8s}")
    print(f"  {'-'*38}")
    for key in ["Hex", "Pen", "Non", "Oct", "Hept"]:
        old = oldCounter.get(key, 0)
        new = newCounter.get(key, 0)
        delta = new - old
        sign = "+" if delta > 0 else ""
        if old > 0 or new > 0:
            print(f"  {key:<12s} {old:>8,} {new:>8,} {sign}{delta:>7,}")

    # 高频预测项
    print(f"\n  [NEW TAGS] Top inferred sugar assignments:")
    allNewTokens = []
    for seq in df["Sugar_Sequence"].dropna():
        tokens = re.findall(
            r'[\w-]+_(?:iupac_parsed|biological_inferred|2d_inferred)',
            str(seq))
        allNewTokens.extend(tokens)
    tagCounter = Counter(allNewTokens)
    for k, v in tagCounter.most_common(15):
        print(f"    {k:40s} {v:>6,}")

    # 剩余死角分析
    stillGeneric = df["Sugar_Sequence"].str.contains(
        genericPat, na=False, regex=True).sum()
    stillTokens = sum(newCounter.values())
    noSmiles = df.loc[
        df["Sugar_Sequence"].str.contains(genericPat, na=False, regex=True)
        & (~df["canonical_smiles"].notna() | (df["canonical_smiles"] == ""))
    ].shape[0]
    noIupac = df.loc[
        df["Sugar_Sequence"].str.contains(genericPat, na=False, regex=True)
        & (~df["iupac_name"].notna() | (df["iupac_name"] == "nan"))
    ].shape[0]

    print(f"\n  [DEAD ENDS] Remaining unresolvable:")
    print(f"    Rows still with generic: {stillGeneric:,}")
    print(f"    Generic tokens remaining: {stillTokens:,}")
    print(f"    No SMILES: {noSmiles:,}")
    print(f"    No IUPAC name: {noIupac:,}")

    # 保存
    df.to_csv(inputPath, index=False, encoding="utf-8-sig")
    print(f"\n  Updated: {inputPath}")
    print(f"  Time: {time.time()-t0:.0f}s")
    print(f"  {'='*70}")


if __name__ == "__main__":
    main()
