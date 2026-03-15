"""
Bio-Value Score — 糖链生物活性价值评分模型
Bio-Value Score — Sugar Chain Bioactivity Value Scoring Model

核心逻辑 (Core Logic):
  对每个 (Sugar_Sequence, Modification) 组合:
    BioValue = Σ pChEMBL_i × target_diversity_weight

  高 pChEMBL 且高靶点多样性 = 高 Bio-Value

依赖: 需先运行 integrate_bioactivity.py 生成 Discovery DB
Dependencies: requires GlycoNP_Final_Discovery_DB.csv

使用方法 (Usage):
  python scripts/compute_bio_value.py [--input PATH] [--mock]
"""
import argparse
import os
import sys

import pandas as pd
import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))


def computeBioValueScores(df: pd.DataFrame) -> pd.DataFrame:
    """
    计算每个 (Sugar_Sequence, Modification_Fingerprint) 的 Bio-Value Score。
    Compute Bio-Value Score for each (Sugar_Sequence, Mod_FP) combination.

    评分模型 (Scoring Model):
      1. potency_score = mean(pChEMBL) — 越高表示越强效
      2. target_diversity = n_unique_targets — 越高表示广谱
      3. frequency = n_compound_activity_pairs — 统计显著性
      4. bio_value = potency_score × log2(1 + target_diversity) × log2(1 + frequency)

    生物学解读: 如果某糖链-修饰组合在多个靶点上都有强活性,
    说明它可能是进化选择的"生物活性增强器"。
    """
    # 确保必要列存在 (Ensure required columns)
    requiredCols = ["Sugar_Sequence", "Glycan_Modifications",
                    "pchembl_value", "target_pref_name"]
    for col in requiredCols:
        if col not in df.columns:
            print(f"  [WARN] Missing column: {col}")
            return pd.DataFrame()

    # 清理 pChEMBL (Clean pChEMBL values)
    df = df.copy()
    df["pchembl_value"] = pd.to_numeric(df["pchembl_value"], errors="coerce")
    df = df[df["pchembl_value"].notna()].copy()

    # 修饰指纹 (Modification fingerprint)
    df["_mod_fp"] = df["Glycan_Modifications"].fillna("None").astype(str)
    df.loc[df["_mod_fp"].isin(["nan", ""]), "_mod_fp"] = "Unmodified"

    # 按 (Sugar_Sequence, _mod_fp) 分组
    grouped = df.groupby(["Sugar_Sequence", "_mod_fp"])

    results = []
    for (sugar, modFp), grp in grouped:
        if not sugar or str(sugar) == "nan":
            continue

        nPairs = len(grp)
        meanPchembl = grp["pchembl_value"].mean()
        maxPchembl = grp["pchembl_value"].max()
        nTargets = grp["target_pref_name"].nunique()
        topTarget = grp["target_pref_name"].value_counts().index[0]

        # Bio-Value Score formula
        bioValue = (meanPchembl
                    * np.log2(1 + nTargets)
                    * np.log2(1 + nPairs))

        results.append({
            "Sugar_Sequence": sugar,
            "Modification": modFp,
            "Bio_Value_Score": round(bioValue, 2),
            "Mean_pChEMBL": round(meanPchembl, 2),
            "Max_pChEMBL": round(maxPchembl, 2),
            "N_Targets": nTargets,
            "N_Pairs": nPairs,
            "Top_Target": topTarget,
        })

    scoreDf = pd.DataFrame(results)
    scoreDf = scoreDf.sort_values("Bio_Value_Score", ascending=False)
    return scoreDf


def main():
    parser = argparse.ArgumentParser(
        description="Bio-Value Score Computation")
    parser.add_argument("--input", type=str, default=None,
                        help="Discovery DB CSV")
    parser.add_argument("--mock", action="store_true",
                        help="Use mock Discovery DB")
    args = parser.parse_args()

    baseDir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    reportDir = os.path.join(baseDir, "reports")

    if args.input:
        inputPath = args.input
    else:
        inputPath = os.path.join(reportDir, "GlycoNP_Final_Discovery_DB.csv")

    print("=" * 60)
    print("  Bio-Value Score Computation")
    print("  糖链生物活性价值评分")
    print("=" * 60)

    if not os.path.exists(inputPath):
        print(f"  [ERROR] Discovery DB not found: {inputPath}")
        print(f"  Run integrate_bioactivity.py first!")
        return

    df = pd.read_csv(inputPath, low_memory=False, dtype=str,
                     encoding="utf-8-sig")
    print(f"  Loaded: {len(df):,} compound-activity pairs")

    scoreDf = computeBioValueScores(df)

    if scoreDf.empty:
        print("  [WARN] No scores computed")
        return

    # 输出报告 (Output report)
    outCsv = os.path.join(reportDir, "Bio_Value_Scores.csv")
    scoreDf.to_csv(outCsv, index=False, encoding="utf-8-sig")
    print(f"\n  Scores saved: {outCsv}")
    print(f"  Total scored combinations: {len(scoreDf)}")

    # Top 10 打印 (Print top 10)
    print(f"\n{'='*60}")
    print(f"  Top 10 Bio-Value Scored Sugar-Modification Combinations")
    print(f"{'='*60}")
    print(f"  {'Sugar':<30} {'Mod':<20} {'Score':>8} {'pCh':>5} "
          f"{'Tgt':>4} {'N':>4} {'Top Target'}")
    print(f"  {'-'*30} {'-'*20} {'-'*8} {'-'*5} {'-'*4} {'-'*4} {'-'*25}")

    for _, row in scoreDf.head(10).iterrows():
        print(f"  {str(row['Sugar_Sequence'])[:30]:<30} "
              f"{str(row['Modification'])[:20]:<20} "
              f"{row['Bio_Value_Score']:>8.1f} "
              f"{row['Mean_pChEMBL']:>5.1f} "
              f"{row['N_Targets']:>4d} "
              f"{row['N_Pairs']:>4d} "
              f"{str(row['Top_Target'])[:25]}")

    # Markdown 报告 (Markdown report)
    mdLines = ["# Bio-Value Score Report",
               "# 糖链生物活性价值评分报告\n",
               f"> Total scored combinations: **{len(scoreDf)}**\n",
               "## Top 20 Sugar–Modification Combinations by Bio-Value\n",
               "| Rank | Sugar Sequence | Modification | Bio-Value | "
               "Mean pChEMBL | Targets | Pairs | Top Target |",
               "|:----:|:--------------|:------------|----------:|"
               "------------:|--------:|------:|:-----------|"]

    for i, (_, row) in enumerate(scoreDf.head(20).iterrows(), 1):
        mdLines.append(
            f"| {i} | `{row['Sugar_Sequence'][:30]}` | "
            f"{row['Modification'][:20]} | "
            f"**{row['Bio_Value_Score']:.1f}** | "
            f"{row['Mean_pChEMBL']:.1f} | "
            f"{row['N_Targets']} | {row['N_Pairs']} | "
            f"{str(row['Top_Target'])[:30]} |"
        )

    mdLines.append(
        "\n---\n> Generated by `compute_bio_value.py` | GlycoNP Pipeline")

    mdPath = os.path.join(reportDir, "Bio_Value_Report.md")
    with open(mdPath, "w", encoding="utf-8") as f:
        f.write("\n".join(mdLines))
    print(f"\n  Report saved: {mdPath}")
    print(f"\n{'='*60}")


if __name__ == "__main__":
    main()
