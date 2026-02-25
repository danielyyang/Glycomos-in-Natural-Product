import os
import sys
import pandas as pd
from rdkit import Chem

# 确保可以引入 lib 文件夹 (Ensure lib folder can be imported)
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from lib.visualizer import StructureVisualizer

def generate_test_images():
    output_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "images", "Test"))
    os.makedirs(output_dir, exist_ok=True)
    
    input_file = os.path.join(os.path.dirname(__file__), "..", "data", "processed", "Coconut_Sugar_Sort.csv")
    df = pd.read_csv(input_file, nrows=100)
    
    viz = StructureVisualizer()
    
    # 手动添加几个合成的测试案例以确保颜色逻辑触发 (Add synthetic test cases for coloring logic)
    test_cases = [
        # O-Glycosylation on Serine (Contains Amino Acid Backbone)
        ("C(C(C(=O)O)N)O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O", "Synthetic_GlycoSerine"),
        # Nucleoside Test (Uridine)
        ("O=C1C=CN([C@H]2O[C@@H](CO)[C@H](O)[C@@H]2O)C(=O)N1", "Synthetic_Uridine")
    ]
    
    print(f"将测试图片输出至 (Outputting test images to): {output_dir}")
    print("-" * 40)
    
    for smiles, name in test_cases:
        out_path = os.path.join(output_dir, f"{name}.png")
        print(f"正在生成 (Generating): {name}.png")
        viz.analyze_glycolipid(smiles, out_path)
    
    # 找几个真实的含氮化合物生成图片 (Find real N-containing compounds to test)
    count = 0
    for idx, row in df.iterrows():
        smiles = row.get('canonical_smiles')
        if pd.notna(smiles) and "N" in str(smiles):
            count += 1
            if count > 5:  # 只生成5个实际样本 (Only generate 5 real samples)
                break
            comp_id = row.get('identifier', f"Row_{idx}")
            out_path = os.path.join(output_dir, f"Test_{comp_id}.png")
            print(f"正在生成 (Generating): Test_{comp_id}.png")
            viz.analyze_glycolipid(smiles, out_path)
            
    print("-" * 40)
    print("✅ 测试图像生成完毕 (Test image generation complete).")

if __name__ == "__main__":
    generate_test_images()
