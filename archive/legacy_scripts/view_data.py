import pandas as pd
import os
import sys

# Fix Windows console emoji encoding
sys.stdout.reconfigure(encoding='utf-8')

def inspect_parquet(file_path):
    """快速读取 Parquet 文件，并打印核心列供终端预览"""
    print(f"🔍 正在读取 Parquet 数据集: {file_path}")
    try:
        # Parquet 极其强大，可以直接秒速加载全量数据而不会轻易 OOM
        df = pd.read_parquet(file_path, engine='pyarrow')
        print(f"\n✅ 数据集加载成功！总行数: {len(df)}")

        # 筛选出我们最关心的几列
        core_columns = ['identifier', 'Glycan_SMILES', 'Aglycan_SMILES', 'Nucleic_Motifs', 'AminoAcid_SMILES', 'Sugar_Sequence', 'murcko_framework']
        existing_cols = [col for col in core_columns if col in df.columns]
        
        print("\n--- 📊 前 5 行核心数据预览 ---")
        print(df[existing_cols].head(5).to_string())
        return df
    except Exception as e:
        print(f"❌ 读取失败: {e}")
        return None

def export_to_csv(df, output_csv, limit=None):
    """将 DataFrame 导出为人类可读的 CSV 文件"""
    if df is None:
        return

    if limit:
        export_df = df.head(limit)
        print(f"\n📦 正在导出前 {limit} 行到 CSV...")
    else:
        export_df = df
        print(f"\n📦 正在导出全量数据到 CSV (警告: 可能会导致大文件)...")
        
    export_df.to_csv(output_csv, index=False, encoding='utf-8-sig')
    print(f"✅ 导出成功！你可以使用 Excel 或文本编辑器打开: {output_csv}")

if __name__ == "__main__":
    base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    # 替换为你的真实 Parquet 文件路径
    parquet_file = os.path.join(base_dir, "reports", "mini_test_output.parquet") # 或者你的全量输出文件

    # 1. 在终端预览数据
    loaded_df = inspect_parquet(parquet_file)
    
    # 2. 为科学家导出一份人类可读的 CSV (比如先导出前 10000 行看看)
    csv_output_path = os.path.join(base_dir, "reports", "human_readable_output_mini.csv")
    export_to_csv(loaded_df, csv_output_path, limit=10000)
