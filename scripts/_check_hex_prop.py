import pandas as pd
import re
from collections import Counter
import os

csv_path = "reports/GlycoNP_Deep_Enriched_Final.csv"
if not os.path.exists(csv_path):
    print(f"Error: {csv_path} not found.")
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

print("\n" + "="*50)
print(f"  糖类标签占比统计 (总独立糖基数: {total:,})")
print("="*50)

# Calculate for the specific generics requested
generics = ["Hex", "Pen", "Hept", "dHex", "HexA"]
for g in generics:
    count = tc.get(g, 0)
    pct = count / total * 100 if total > 0 else 0
    print(f"  {g:<15} {count:>8,}   ({pct:5.2f}%)")

print("\n  Top 10 具体糖 (Specific Sugars):")
print("-" * 50)
specific_count = 0
rank = 1
for sugar, count in tc.most_common():
    if sugar not in ["Hex", "Pen", "Hept", "dHex", "HexA", "Non", "Oct", "HexN", "HexNAc"]:
        pct = count / total * 100 if total > 0 else 0
        if rank <= 10:
            print(f"  {rank:2d}. {sugar:<25} {count:>8,} ({pct:5.2f}%)")
        specific_count += count
        rank += 1

specific_pct = specific_count / total * 100 if total > 0 else 0
generic_count = total - specific_count
generic_pct = generic_count / total * 100 if total > 0 else 0

print("-" * 50)
print(f"  总泛指标签 (Generics): {generic_count:>8,} ({generic_pct:5.2f}%)")
print(f"  总具体糖 (Specifics): {specific_count:>8,} ({specific_pct:5.2f}%)")
print("="*50)
