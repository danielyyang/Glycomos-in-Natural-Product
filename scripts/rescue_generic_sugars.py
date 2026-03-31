"""
GlycoNP 泛指糖精确修复管线 (Generic Sugar Rescue Pipeline)

将 Hex/Pen/dHex/Hept/Oct/Non 等泛化标签修复为精确单糖名称。
Uses ONLY 1 strategy:
  A. Text-Mining from compound name/iupac_name/synonyms
  (Strategies C and D based on statistical scaffold prior have been stripped 
  to ensure strict chemical compliance and avoid false positive hallucinations.)

使用方法 (Usage):
  python scripts/rescue_generic_sugars.py
"""
import os
import re
import sys
import time
from collections import Counter
from typing import Optional

import pandas as pd

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

# =====================================================================
# 策略 A: 化合物名称文本挖掘字典 (Text-Mining Dictionary)
# 从化合物名称中提取明确的糖名称线索
# =====================================================================
# 优先级排序: 最长匹配优先, 避免 "glucoside" 误匹配 "galactoside"
TEXT_MINING_RULES = [
    # 精确多字匹配 (Multi-word exact, highest priority)
    (r'glucuronic\s*acid|glucuronide|glucuronosyl', 'D-GlcA'),
    (r'galacturonic\s*acid|galacturonide', 'D-GalA'),
    (r'N-acetylglucosamin|GlcNAc', 'D-GlcNAc'),
    (r'N-acetylgalactosamin|GalNAc', 'D-GalNAc'),
    (r'neuraminic|sialic|Neu5Ac|NeuAc', 'Neu5Ac'),

    # 精确单糖名 + 复合词根匹配 (Exact names + compound word roots)
    (r'glucopyranosid|glucosinolat|glucosid|glucosyl|glucofuranos|gluco(?:se)?(?![a-z])', 'D-Glc'),
    (r'galactopyranosid|galactolipid|galactosid|galactosyl|galacto(?:se)?(?![a-z])', 'D-Gal'),
    (r'mannopyranosid|mannoprotein|mannosid|mannosyl|manno(?:se)?(?![a-z])', 'D-Man'),
    (r'xylopyranosid|xylosid|xylosyl|xylan|xylo(?:se)?(?![a-z])', 'D-Xyl'),
    (r'arabinopyranosid|arabinosid|arabinosyl|arabinan|arabino(?:se)?(?![a-z])', 'L-Ara'),
    (r'rhamnopyranosid|rhamnosid|rhamnosyl|rhamno(?:se)?(?![a-z])', 'L-Rha'),
    (r'fucopyranosid|fucosid|fucosyl|fucoidan|fuco(?:se)?(?![a-z])', 'L-Fuc'),
    (r'ribopyranosid|ribosid|ribosyl|ribo(?:se)?(?![a-z])', 'D-Rib'),
    (r'fructopyranosid|fructosid|fructosyl|fructan|fructo(?:se)?(?![a-z])', 'D-Fru'),

    # 脱氧糖特异名 (Deoxy sugar specific names)
    (r'digitalosid|digitalose', 'D-Dig'),
    (r'oleandrosid|oleandrose', 'L-Ole'),
    (r'cymarosid|cymarose', 'D-Cym'),
    (r'thevetosid|thevetose', 'D-The'),
    (r'boivinosid|boivinose', 'D-Boi'),
    (r'quinovo(?:se|sid)', 'D-Qui'),
]
# 预编译正则 (Pre-compile regex)
COMPILED_TEXT_RULES = [
    (re.compile(pattern, re.IGNORECASE), sugar)
    for pattern, sugar in TEXT_MINING_RULES
]


def strategyA_textMining(
    name: Optional[str],
    iupacName: Optional[str],
    synonyms: Optional[str],
) -> Optional[str]:
    """
    策略 A: 从化合物名称/IUPAC名/同义词中文本挖掘精确糖名。
    Strategy A: Text-mine exact sugar name from compound name fields.

    Returns:
        精确糖名 or None (if no match found)
    """
    # 合并所有可用文本线索
    textPool = " ".join(filter(None, [
        str(name) if name and str(name) != "nan" else None,
        str(iupacName) if iupacName and str(iupacName) != "nan" else None,
        str(synonyms) if synonyms and str(synonyms) != "nan" else None,
    ]))

    if not textPool.strip():
        return None

    for regex, sugar in COMPILED_TEXT_RULES:
        if regex.search(textPool):
            return f"{sugar}_predicted" 

    return None


# =====================================================================
# 主管线: 逐行修复 (Main Pipeline)
# =====================================================================
def rescueSingleToken(
    token: str,
    name: Optional[str],
    iupacName: Optional[str],
    synonyms: Optional[str],
    logCounter: dict,
) -> str:
    """
    对单个泛化 token 执行修复, 仅使用严谨的文本挖掘 (Strategy A)。
    Rescue single generic token via strict Text-Mining.
    """
    # 不是泛化标签, 直接返回
    if token not in ("Hex", "Pen", "dHex", "HexA", "Non", "Oct", "Hept"):
        return token

    # 策略 A: 文本匹配
    resultA = strategyA_textMining(name, iupacName, synonyms)
    if resultA:
        logCounter["A"] += 1
        return resultA

    # 全部失败, 保留原标签
    logCounter["MISS"] += 1
    return token


def rescueSequence(
    seq: str,
    name: Optional[str],
    iupacName: Optional[str],
    synonyms: Optional[str],
    logCounter: dict,
) -> str:
    """
    修复整条 Sugar_Sequence 中的泛化 token。
    Rescue all generic tokens in a Sugar_Sequence string.
    """
    genericPattern = re.compile(r'\b(Hex|Pen|dHex|HexA|Non|Oct|Hept)\b')

    def replaceMatch(m):
        token = m.group(0)
        return rescueSingleToken(
            token, name, iupacName, synonyms, logCounter)

    return genericPattern.sub(replaceMatch, seq)


def main():
    baseDir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    reportDir = os.path.join(baseDir, "reports")
    inputPath = os.path.join(reportDir, "GlycoNP_Pipeline_Full_Cleaned.csv")

    print("=" * 70)
    print("  GlycoNP Generic Sugar Rescue Pipeline (STRICT NLP ONLY)")
    print("  泛指糖精确修复管线 (纯文本挖掘模式)")
    print("=" * 70)

    t0 = time.time()
    if not os.path.exists(inputPath):
        print(f"Error: {inputPath} not found.")
        sys.exit(1)

    df = pd.read_csv(inputPath, low_memory=False, encoding="utf-8-sig")
    total = len(df)
    print(f"  Loaded: {total:,} rows ({time.time()-t0:.1f}s)")

    # ===== Step 1: 审计 =====
    genericPattern = r'\bHex\b|\bPen\b|\bdHex\b|\bHexA\b|\bNon\b|\bOct\b|\bHept\b'
    targetMask = df["Sugar_Sequence"].str.contains(
        genericPattern, na=False, regex=True)
    targetCount = targetMask.sum()
    print(f"\n  [AUDIT] Rows with generic sugars: {targetCount:,}")

    # 旧分布 (Old distribution)
    oldTokens = []
    for seq in df.loc[targetMask, "Sugar_Sequence"].dropna():
        oldTokens.extend(re.findall(genericPattern, str(seq)))
    oldCounter = Counter(oldTokens)
    print(f"  Generic tokens breakdown:")
    for k, v in oldCounter.most_common():
        print(f"    {k:10s} {v:>6,}")

    # ===== Step 2: 执行修复 =====
    print(f"\n  [EXECUTE] Rescuing {targetCount:,} rows using Text-Mining ONLY...")
    logCounter = {"A": 0, "MISS": 0}

    df["Sugar_Sequence_Before"] = df["Sugar_Sequence"].copy()

    rescued = 0
    for idx in df.index[targetMask]:
        oldSeq = str(df.at[idx, "Sugar_Sequence"])
        newSeq = rescueSequence(
            oldSeq,
            df.at[idx, "name"] if pd.notna(df.at[idx, "name"]) else None,
            df.at[idx, "iupac_name"] if pd.notna(df.at[idx, "iupac_name"]) else None,
            df.at[idx, "synonyms"] if pd.notna(df.at[idx, "synonyms"]) else None,
            logCounter,
        )
        if newSeq != oldSeq:
            df.at[idx, "Sugar_Sequence"] = newSeq
            rescued += 1

    print(f"\n  [RESULT] Rescue Summary:")
    print(f"    Total rows modified: {rescued:,} / {targetCount:,}")
    print(f"    Strategy A (Text-Mining):   {logCounter['A']:>6,} tokens")
    print(f"    Unrescued (MISS):           {logCounter['MISS']:>6,} tokens")

    totalTokensRescued = logCounter["A"]
    totalTokens = totalTokensRescued + logCounter["MISS"]
    rescueRate = totalTokensRescued / totalTokens * 100 if totalTokens > 0 else 0
    print(f"    Rescue rate: {rescueRate:.1f}%")

    # ===== 新分布 =====
    newTokens = []
    for seq in df["Sugar_Sequence"].dropna():
        newTokens.extend(re.findall(genericPattern, str(seq)))
    newCounter = Counter(newTokens)

    print(f"\n  [COMPARISON] Before → After:")
    print(f"  {'Token':<12s} {'Before':>8s} {'After':>8s} {'Delta':>8s}")
    print(f"  {'-'*38}")
    for key in ["Hex", "Pen", "dHex", "HexA", "Non", "Oct", "Hept"]:
        old = oldCounter.get(key, 0)
        new = newCounter.get(key, 0)
        delta = new - old
        if old > 0 or new > 0:
            sign = "+" if delta > 0 else ""
            print(f"  {key:<12s} {old:>8,} {new:>8,} {sign}{delta:>7,}")

    # 新增精确糖统计 (Newly assigned precise sugars)
    print(f"\n  [NEW ASSIGNMENTS] Top precise sugars assigned:")
    allNewTokens = []
    for seq in df.loc[targetMask, "Sugar_Sequence"].dropna():
        tokens = re.findall(
            r'Neu5Ac|Neu5Gc|KDO|'
            r'[DL]-[A-Z][a-z]+[A-Z]?[a-z]*(?:_predicted)?',
            str(seq))
        allNewTokens.extend(tokens)
    assignedCounter = Counter(allNewTokens)
    for k, v in assignedCounter.most_common(15):
        print(f"    {k:25s} {v:>6,}")

    # ===== 保存 =====
    df.drop(columns=["Sugar_Sequence_Before"], inplace=True)
    outPath = os.path.join(reportDir, "GlycoNP_Pipeline_Full_Cleaned.csv")
    df.to_csv(outPath, index=False, encoding="utf-8-sig")
    print(f"\n  Updated: {outPath}")
    print(f"  Time: {time.time()-t0:.0f}s")
    print(f"  {'='*70}")


if __name__ == "__main__":
    main()
