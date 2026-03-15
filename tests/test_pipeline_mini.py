import os
import sys
import pandas as pd

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from scripts.batch_processor import process_chunk

def run_mini_batch_test():
    base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    input_csv = os.path.join(base_dir, "data", "Coconut.csv")
    output_parquet = os.path.join(base_dir, "reports", "mini_test_output.parquet")
    
    print("--- [START] Starting Mini-Batch Dry Run (First 1000 Rows) ---")
    
    # 1. 抽取微缩样本 (Mini-Batch)
    try:
        df_raw = pd.read_csv(input_csv, nrows=1000, dtype=str, encoding='utf-8-sig', low_memory=False)
        print(f"Loaded {len(df_raw)} rows from {input_csv}.")
    except Exception as e:
        print(f"Failed to load raw data: {e}")
        return

    # 2. 贯通 7-Phase 引擎群
    processed_df = process_chunk(df_raw)
    
    if processed_df.empty:
        print("No sugar-containing molecules found in the first 1000 rows.")
        return
        
    print(f"Pipeline processed successfully. {len(processed_df)} sugar-containing molecules identified.")
    
    # 确保列名准确 (应对 Phase 3 核酸字段名的微小变更)
    if 'NUCLEOTIDES_SMILES' in processed_df.columns:
        processed_df.rename(columns={'NUCLEOTIDES_SMILES': 'Nucleic_Motifs'}, inplace=True)
    
    # 3. 极简落盘
    import pyarrow as pa
    import pyarrow.parquet as pq
    
    # 确保全为字符串类型
    processed_df = processed_df.astype(str)
    table = pa.Table.from_pandas(processed_df)
    pq.write_table(table, output_parquet)
    print(f"Saved mini-batch to {output_parquet}")
    
    # 4. 数据探伤 (Data Profiling)
    print("\n--- [SEARCH] Data Profiling (Loading from Parquet) ---")
    df_check = pd.read_parquet(output_parquet)
    
    # 根据用户要求的关键列进行过滤打印
    target_columns = ['identifier', 'Glycan_SMILES', 'Aglycan_SMILES', 'Nucleic_Motifs', 'AminoAcid_SMILES', 'Sugar_Sequence', 'murcko_framework']
    
    # 确保列都存在，不存在赋 NULL (防意外)
    for col in target_columns:
        if col not in df_check.columns:
            df_check[col] = "MISSING_COL"
            
    print(df_check[target_columns].head(5).to_string(index=False))
    print("\n[SUCCESS] Mini-batch dry run completed.")

if __name__ == "__main__":
    run_mini_batch_test()
