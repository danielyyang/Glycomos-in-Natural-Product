"""
GlycoNP 糖序列重新匹配 — 使用升级后的库 (Sugar Sequence Re-Matching)

在全量清洗数据上重新运行 sugar_sequence.analyze_glycan(),
使用更新的精确库 (Neu5Ac, KDO, Heptose) + 碳计数退避 + 开环过滤。

使用方法 (Usage):
  python scripts/rematch_sugar_sequence.py
"""
import os
import sys
import time
import re
from collections import Counter

import pandas as pd
from tqdm import tqdm
from rdkit import Chem

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "lib")))

from lib import sugar_sequence


def main():
    baseDir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    reportDir = os.path.join(baseDir, "reports")
    inputPath = os.path.join(reportDir, "GlycoNP_Pipeline_Full_Cleaned.csv")

    print("=" * 70)
    print("  Sugar Sequence Re-Matching (Upgraded Library)")
    print("  糖序列重新匹配 (升级库版)")
    print("=" * 70)
    print(f"  Input: {inputPath}")

    t0 = time.time()
    df = pd.read_csv(inputPath, low_memory=False, encoding="utf-8-sig")
    total = len(df)
    print(f"  Loaded: {total:,} rows ({time.time()-t0:.1f}s)")

    # 保存旧序列用于对比 (Save old sequences for comparison)
    df["Sugar_Sequence_Old"] = df["Sugar_Sequence"].copy()

    # 收集旧统计 (Old statistics)
    oldTokens = []
    for seq in df["Sugar_Sequence_Old"].dropna().astype(str):
        oldTokens.extend(re.findall(
            r'Neu5Ac|Neu5Gc|KDO|Hept|Oct|Non|Non_Cyclic_Invalid|'
            r'[DL]-[A-Z][a-z]+[A-Z]?[a-z]*|Hex|dHex|Pen|HexA|HexN|'
            r'HexNAc|Unk|Invalid|Error|Unknown',
            seq))
    oldCounter = Counter(oldTokens)
    print(f"\n  Old Sugar_Sequence stats:")
    print(f"    Total tokens: {len(oldTokens):,}")
    print(f"    Hex: {oldCounter.get('Hex', 0):,}")
    print(f"    Unknown: {oldCounter.get('Unknown', 0):,}")

    # 重新匹配 (Re-match)
    print(f"\n  Re-matching {total:,} molecules...")
    rematched = 0
    errors = 0

    for idx in tqdm(df.index, desc="  Re-matching", ncols=80):
        smiles = str(df.at[idx, "canonical_smiles"])
        if smiles in ("nan", "", "None"):
            continue

        try:
            newSeq, newMods = sugar_sequence.analyze_glycan(smiles)
            if newSeq:
                df.at[idx, "Sugar_Sequence"] = newSeq
                if newMods:
                    df.at[idx, "Glycan_Modifications"] = newMods
                rematched += 1
        except Exception:
            errors += 1

    print(f"\n  Re-matched: {rematched:,} / {total:,}")
    print(f"  Errors: {errors}")

    # ================================================================
    # Stereo-Rescue Pass (立体化学确证 — 仅对含泛指标签的行执行)
    # 结果写入新列 Rescued_Sugar_Sequence, 不覆盖原 Sugar_Sequence
    # Results go to new column, preserving original Sugar_Sequence
    # ================================================================
    print(f"\n  === Stereo-Rescue Pass ===")
    try:
        from lib.stereo_rescue import rescueSugarSequence, GENERIC_SUGAR_LABELS
        rescueCount = 0
        rescueNameCount = 0
        rescueChemblCount = 0

        # 初始化新列 (Initialize new columns)
        df["Rescued_Sugar_Sequence"] = ""
        df["Rescue_Method"] = ""

        # 找出含泛指标签的行 (Find rows with generic sugar labels)
        genericMask = df["Sugar_Sequence"].astype(str).apply(
            lambda s: any(label in s for label in GENERIC_SUGAR_LABELS)
        )
        genericIndices = df.index[genericMask]
        print(f"  Rows with generic labels: {len(genericIndices):,}")

        for idx in tqdm(genericIndices, desc="  Stereo-Rescue", ncols=80):
            smiles = str(df.at[idx, "canonical_smiles"])
            currSeq = str(df.at[idx, "Sugar_Sequence"])
            name = str(df.at[idx, "name"]) if "name" in df.columns else None
            iupacName = str(df.at[idx, "iupac_name"]) if "iupac_name" in df.columns else None

            rescuedSeq, method = rescueSugarSequence(
                smiles, currSeq, name, iupacName
            )

            if method and rescuedSeq != currSeq:
                df.at[idx, "Rescued_Sugar_Sequence"] = rescuedSeq
                df.at[idx, "Rescue_Method"] = method
                rescueCount += 1
                if "Name" in method:
                    rescueNameCount += 1
                elif "ChEMBL" in method:
                    rescueChemblCount += 1

        print(f"  Rescued: {rescueCount:,} / {len(genericIndices):,}")
        print(f"    Name-Inferred: {rescueNameCount:,}")
        print(f"    ChEMBL-Candidate: {rescueChemblCount:,}")
    except Exception as e:
        print(f"  Stereo-Rescue skipped (error: {e})")

    # 新统计 (New statistics)
    newTokens = []
    for seq in df["Sugar_Sequence"].dropna().astype(str):
        newTokens.extend(re.findall(
            r'Neu5Ac|Neu5Gc|KDO|Hept|Oct|Non|Non_Cyclic_Invalid|'
            r'[DL]-[A-Z][a-z]+[A-Z]?[a-z]*|Hex|dHex|Pen|HexA|HexN|'
            r'HexNAc|Unk|Invalid|Error|Unknown',
            seq))
    newCounter = Counter(newTokens)

    print(f"\n  {'='*60}")
    print(f"  Comparison Report (Old → New)")
    print(f"  {'='*60}")
    print(f"  Total tokens: {len(oldTokens):,} → {len(newTokens):,}")

    # 关键变化 (Key changes)
    compareKeys = [
        "Hex", "dHex", "Pen", "Unknown", "Invalid", "Error",
        "Hept", "Oct", "Non", "Non_Cyclic_Invalid",
        "Neu5Ac", "Neu5Gc", "KDO", "L-D-Hep", "D-D-Hep",
        "D-Glc", "D-Gal", "D-Man", "L-Rha", "L-Fuc",
        "D-Xyl", "L-Ara", "D-GlcA", "D-GalA",
    ]

    print(f"\n  {'Token':<25s} {'Old':>8s} {'New':>8s} {'Delta':>8s}")
    print(f"  {'-'*51}")
    for key in compareKeys:
        old = oldCounter.get(key, 0)
        new = newCounter.get(key, 0)
        delta = new - old
        if old > 0 or new > 0:
            sign = "+" if delta > 0 else ""
            print(f"  {key:<25s} {old:>8,} {new:>8,} {sign}{delta:>7,}")

    # 保存 (Save)
    df.drop(columns=["Sugar_Sequence_Old"], inplace=True)
    outPath = os.path.join(reportDir, "GlycoNP_Pipeline_Full_Cleaned.csv")
    df.to_csv(outPath, index=False, encoding="utf-8-sig")
    print(f"\n  Updated: {outPath}")
    print(f"  Time: {time.time()-t0:.0f}s")
    print(f"  {'='*60}")


if __name__ == "__main__":
    main()

