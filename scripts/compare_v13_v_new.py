import pandas as pd
import os

print("Loading V13...")
try:
    df_v13 = pd.read_csv("reports/GlycoNP_Deep_Enriched_v13_Final.csv", low_memory=False)
except FileNotFoundError:
    print("V13 not found at reports/GlycoNP_Deep_Enriched_v13_Final.csv")
    exit()

print("Loading New Final...")
try:
    df_new = pd.read_csv("reports/GlycoNP_Deep_Enriched_Final.csv", low_memory=False)
except FileNotFoundError:
    print("New Final not found at reports/GlycoNP_Deep_Enriched_Final.csv")
    exit()

# Merge on identifier
df = pd.merge(df_v13[["identifier", "canonical_smiles", "name", "Sugar_Sequence", "Consensus_Sugar_Sequence"]], 
              df_new[["identifier", "Sugar_Sequence", "Consensus_Sugar_Sequence"]], 
              on="identifier", 
              suffixes=('_v13', '_new'))

# Find cases where v13 had specific sugar (e.g. D-Glc) and new has Hex
# We use Consensus_Sugar_Sequence for _new to see what the raw engine output prior to anti-regression backfill
mask = (df["Sugar_Sequence_v13"].str.contains("D-Glc", na=False)) & \
       (df["Sugar_Sequence_new"].str.contains("Hex", na=False))

diff = df[mask].head(10)

print(f"Found {len(df[mask])} rows where V13 was D-Glc but New is Hex.")
print("="*80)
for idx, row in diff.iterrows():
    print(f"ID:      {row['identifier']}")
    print(f"Name:    {row['name']}")
    print(f"SMILES:  {row['canonical_smiles']}")
    print(f"V13 Seq: {row['Sugar_Sequence_v13']}")
    print(f"New Seq: {row['Sugar_Sequence_new']}")
    print(f"Has @@:  {'@@' in str(row['canonical_smiles'])}")
    print("-" * 80)
