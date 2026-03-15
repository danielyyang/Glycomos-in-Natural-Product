"""
Phase 6 Multi-Sheet Excel Export — 按大类分 Sheet 导出
Phase 6 Multi-Sheet Excel Export — Export by Superclass

将 9.4 万条数据按 Superclass 拆分到不同 Excel Sheet。
Split 94K compounds into separate Excel sheets by Superclass.

使用方法 (Usage):
  python scripts/export_multisheet_excel.py [--input PATH]
"""
import argparse
import os
import sys

import pandas as pd

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))


# 需要导出的核心列 (Core columns for export)
EXPORT_COLUMNS = [
    "name", "standard_inchi_key", "molecular_formula", "molecular_weight",
    "canonical_smiles", "Glycan_SMILES", "Aglycan_SMILES",
    "Sugar_Sequence", "Glycan_Modifications",
    "Has_Nucleotide", "Has_Peptide",
    "Murcko_Scaffold", "Superclass",
    "organisms", "organism_taxonomy_05family",
    "alogp", "np_likeness", "dois",
]


def cleanSuperclass(val: str) -> str:
    """清理 Superclass 名称为合法 Sheet 名。"""
    if not val or val == "nan":
        return "Unclassified"
    s = str(val).strip()
    if "(Tanimoto=" in s:
        s = s[:s.index("(Tanimoto=")].strip()
    if s.startswith("Glycolipid"):
        s = "Glycolipids"
    # Excel Sheet 名限制 31 字符, 不能含 [\]/*?:
    s = s.replace("/", "-").replace("\\", "-").replace("*", "")
    s = s.replace("[", "(").replace("]", ")").replace("?", "").replace(":", "-")
    return s[:31] if s else "Unclassified"


def main():
    parser = argparse.ArgumentParser(
        description="Multi-Sheet Excel Export by Superclass")
    parser.add_argument("--input", type=str, default=None, help="Input CSV")
    args = parser.parse_args()

    baseDir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    reportDir = os.path.join(baseDir, "reports")

    if args.input:
        inputPath = args.input
    else:
        inputPath = os.path.join(reportDir, "GlycoNP_Pipeline_Full.csv")

    print("=" * 60)
    print("  Multi-Sheet Excel Export by Superclass")
    print("  按大类分 Sheet 导出")
    print("=" * 60)
    print(f"  Input: {inputPath}")

    df = pd.read_csv(inputPath, low_memory=False, dtype=str,
                     encoding="utf-8-sig")
    total = len(df)
    print(f"  Loaded: {total:,} rows")

    # 清理 Superclass (Clean Superclass names)
    df["_sheet"] = df["Superclass"].apply(cleanSuperclass)

    # 选择可用列 (Select available columns)
    availCols = [c for c in EXPORT_COLUMNS if c in df.columns]
    print(f"  Export columns: {len(availCols)}")

    outputPath = os.path.join(reportDir, "GlycoNP_MultiSheet.xlsx")

    with pd.ExcelWriter(outputPath, engine="openpyxl") as writer:
        # Summary sheet first
        summaryData = df["_sheet"].value_counts().reset_index()
        summaryData.columns = ["Superclass", "Count"]
        summaryData["Percentage"] = (summaryData["Count"] / total * 100).round(1)
        summaryData.to_excel(writer, sheet_name="Summary", index=False)
        print(f"\n  Sheet 'Summary': {len(summaryData)} classes")

        # Per-class sheets
        sheetCount = 0
        for sheetName in summaryData["Superclass"].tolist():
            subset = df[df["_sheet"] == sheetName][availCols].copy()

            # 按 Sugar_Sequence 频率排序 (Sort by sugar frequency)
            if "Sugar_Sequence" in subset.columns:
                freqMap = subset["Sugar_Sequence"].value_counts().to_dict()
                subset["_freq"] = subset["Sugar_Sequence"].map(freqMap).fillna(0)
                subset = subset.sort_values("_freq", ascending=False)
                subset = subset.drop(columns=["_freq"])

            # Excel sheet name max 31 chars
            safeSheetName = sheetName[:31]
            try:
                subset.to_excel(writer, sheet_name=safeSheetName, index=False)
                sheetCount += 1
                print(f"    Sheet '{safeSheetName}': {len(subset):,} rows")
            except Exception as e:
                print(f"    [WARN] Skip '{safeSheetName}': {e}")

    fileSize = os.path.getsize(outputPath) / 1024 / 1024
    print(f"\n  Output: {outputPath}")
    print(f"  File size: {fileSize:.1f} MB")
    print(f"  Total sheets: {sheetCount + 1} (including Summary)")
    print(f"\n{'='*60}")


if __name__ == "__main__":
    main()
