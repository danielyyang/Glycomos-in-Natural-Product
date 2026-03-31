import pandas as pd
import re
from collections import Counter

df = pd.read_csv("reports/GlycoNP_Deep_Enriched_Final.csv", low_memory=False)

seqs = df["Sugar_Sequence"].dropna().astype(str).tolist()

PAT = (r'Neu5Ac|Neu5Gc|KDO|'
       r'[DL]-[A-Z][a-z]+[A-Z]?[a-z]*(?:\([^)]*\))?|'
       r'Hex|dHex|Pen|HexA|HexN|HexNAc|Non|Oct|Hept')
allTokens = []
for seq in seqs:
    allTokens.extend(re.findall(PAT, seq))

tc = Counter(allTokens)
total = sum(tc.values())

for rank, (sugar, count) in enumerate(tc.most_common(20), 1):
    pct = count / total * 100 if total > 0 else 0
    print(f"{rank:2d}. {sugar:<35} {count:>8,} ({pct:5.1f}%)")

genericCount = sum(v for k, v in tc.items() 
                   if k in ("Hex","Pen","dHex","HexA","Non","Oct","Hept"))
print(f"\nGeneric: {genericCount:>8,}")
print(f"Total:   {total:>8,}")
