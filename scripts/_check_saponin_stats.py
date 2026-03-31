import pandas as pd
import re
from collections import Counter
import os

csv_path = "reports/GlycoNP_Saponin_DB.csv"
if not os.path.exists(csv_path):
    print(f"Error: {csv_path} not found. Please ensure extract_saponins.py has finished.")
    exit()

df = pd.read_csv(csv_path, low_memory=False)

if "Sugar_Sequence" not in df.columns:
    print("Error: 'Sugar_Sequence' column not found.")
    exit()

seqs = df["Sugar_Sequence"].dropna().astype(str).tolist()

PAT = (r'Neu5Ac|Neu5Gc|KDO|'
       r'[DL]-[A-Z][a-z]+[A-Z]?[a-z]*(?:\([^)]*\))?|'
       r'Hex|dHex|Pen|HexA|HexN|HexNAc|Non|Oct|Hept')
allTokens = []
for seq in seqs:
    allTokens.extend(re.findall(PAT, seq))

tc = Counter(allTokens)
total = sum(tc.values())

print("\n" + "="*55)
print(f"  🌟 皂苷专属子集 (Saponin Database) 数据洞察 🌟")
print("="*55)
print(f"  总样本数 (Total Saponin Molecules): {len(df):,}")
print(f"  包含的总糖基数 (Total Sugar Units):   {total:,}")
avg_sugar = total / len(df) if len(df) > 0 else 0
print(f"  平均糖链长度 (Avg Sugar per Saponin): {avg_sugar:.2f}糖/分子")
print("-" * 55)

# Calculate for the specific generics requested
generics = ["Hex", "Pen", "Hept", "dHex", "HexA"]
print("\n  [泛指标签分布 (Generics)]:")
for g in generics:
    count = tc.get(g, 0)
    pct = count / total * 100 if total > 0 else 0
    print(f"    {g:<15} {count:>8,}   ({pct:5.2f}%)")

print("\n  [Top 10 专属单糖 (Specific Sugars)]:")
specific_count = 0
rank = 1
for sugar, count in tc.most_common():
    if sugar not in ["Hex", "Pen", "Hept", "dHex", "HexA", "Non", "Oct", "HexN", "HexNAc"]:
        pct = count / total * 100 if total > 0 else 0
        if rank <= 10:
            print(f"    {rank:2d}. {sugar:<25} {count:>8,} ({pct:5.2f}%)")
        specific_count += count
        rank += 1

specific_pct = specific_count / total * 100 if total > 0 else 0
generic_count = total - specific_count
generic_pct = generic_count / total * 100 if total > 0 else 0

print("-" * 55)
print(f"  总泛指标签比例 (Generics Rate):  {generic_count:>8,} ({generic_pct:5.2f}%)")
print(f"  极其纯净的精确糖 (Specifics Rate): {specific_count:>8,} ({specific_pct:5.2f}%)")
print("="*55)
