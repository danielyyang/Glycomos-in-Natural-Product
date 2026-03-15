"""
ChEMBL SQLite → CSV 提取工具
ChEMBL SQLite → CSV Extraction Tool

从 ChEMBL SQLite 数据库中提取 integrate_bioactivity.py 所需的两张表:
  1. chembl_molecules.csv (standard_inchi_key → molecule_chembl_id)
  2. chembl_activities.csv (chembl_id → target, IC50, pChEMBL)

使用方法 (Usage):
  1. 从 https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/
     下载 chembl_36_sqlite.tar.gz (~4GB)
  2. 解压得到 chembl_36.db (SQLite 文件)
  3. 运行:
     python scripts/extract_chembl_tables.py --db path/to/chembl_36.db

输出文件自动保存到 data/chembl/ 目录。
"""
import argparse
import os
import sqlite3
import sys
import time

import pandas as pd

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))


def extractMolecules(conn: sqlite3.Connection, outputPath: str) -> int:
    """
    提取分子映射表: standard_inchi_key → molecule_chembl_id。
    Extract molecule lookup: standard_inchi_key → molecule_chembl_id.

    SQL 逻辑: 从 molecule_dictionary 和 compound_structures 联合查询,
    仅保留有 InChIKey 的分子。
    """
    print("\n  [1/2] Extracting molecule dictionary...")
    t0 = time.time()

    query = """
    SELECT
        cs.standard_inchi_key,
        md.chembl_id AS molecule_chembl_id
    FROM
        molecule_dictionary md
    JOIN
        compound_structures cs ON md.molregno = cs.molregno
    WHERE
        cs.standard_inchi_key IS NOT NULL
        AND cs.standard_inchi_key != ''
    """

    df = pd.read_sql_query(query, conn)
    df = df.drop_duplicates(subset=["standard_inchi_key"])
    df.to_csv(outputPath, index=False, encoding="utf-8-sig")

    elapsed = time.time() - t0
    print(f"    Rows: {len(df):,} ({elapsed:.1f}s)")
    print(f"    Saved: {outputPath}")
    return len(df)


def extractActivities(conn: sqlite3.Connection, outputPath: str) -> int:
    """
    提取活性数据: chembl_id → 靶点/IC50/pChEMBL。
    Extract activities: chembl_id → target/IC50/pChEMBL.

    SQL 逻辑: 多表联查 activities → assays → target_dictionary,
    仅保留有 pChEMBL 值的高质量数据 (Binding assay 优先)。
    过滤条件:
      - standard_type IN ('IC50', 'Ki', 'EC50', 'Kd')
      - pchembl_value IS NOT NULL (保证数据可比较)
    """
    print("\n  [2/2] Extracting activities (this may take a few minutes)...")
    t0 = time.time()

    query = """
    SELECT
        md.chembl_id AS molecule_chembl_id,
        td.pref_name AS target_pref_name,
        td.organism AS target_organism,
        act.standard_type,
        act.standard_value,
        act.standard_units,
        act.pchembl_value,
        ass.assay_type
    FROM
        activities act
    JOIN
        molecule_dictionary md ON act.molregno = md.molregno
    JOIN
        assays ass ON act.assay_id = ass.assay_id
    JOIN
        target_dictionary td ON ass.tid = td.tid
    WHERE
        act.standard_type IN ('IC50', 'Ki', 'EC50', 'Kd')
        AND act.pchembl_value IS NOT NULL
        AND td.pref_name IS NOT NULL
    """

    # 分块写入以控制内存 (Chunked write for memory control)
    chunkSize = 500_000
    totalRows = 0
    firstChunk = True

    for chunk in pd.read_sql_query(query, conn, chunksize=chunkSize):
        if firstChunk:
            chunk.to_csv(outputPath, index=False, encoding="utf-8-sig", mode="w")
            firstChunk = False
        else:
            chunk.to_csv(outputPath, index=False, encoding="utf-8-sig",
                         mode="a", header=False)
        totalRows += len(chunk)
        print(f"    ... {totalRows:,} rows extracted", end="\r")

    elapsed = time.time() - t0
    print(f"    Rows: {totalRows:,} ({elapsed:.1f}s)          ")
    print(f"    Saved: {outputPath}")
    return totalRows


def main():
    parser = argparse.ArgumentParser(
        description="Extract ChEMBL SQLite → CSV for GlycoNP integration")
    parser.add_argument("--db", type=str, required=True,
                        help="Path to ChEMBL SQLite database file")
    parser.add_argument("--output-dir", type=str, default=None,
                        help="Output directory (default: data/chembl/)")
    args = parser.parse_args()

    if not os.path.exists(args.db):
        print(f"  [ERROR] Database not found: {args.db}")
        print(f"  Download from: https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/")
        return

    baseDir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    outputDir = args.output_dir or os.path.join(baseDir, "data", "chembl")
    os.makedirs(outputDir, exist_ok=True)

    molPath = os.path.join(outputDir, "chembl_molecules.csv")
    actPath = os.path.join(outputDir, "chembl_activities.csv")

    print("=" * 60)
    print("  ChEMBL SQLite → CSV Extraction")
    print("=" * 60)
    print(f"  Database: {args.db}")
    print(f"  Output: {outputDir}")

    conn = sqlite3.connect(args.db)

    try:
        nMol = extractMolecules(conn, molPath)
        nAct = extractActivities(conn, actPath)
    finally:
        conn.close()

    print(f"\n{'='*60}")
    print(f"  Extraction Complete!")
    print(f"    Molecules: {nMol:,} → {molPath}")
    print(f"    Activities: {nAct:,} → {actPath}")
    print(f"\n  Next step:")
    print(f"    python scripts/integrate_bioactivity.py")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
