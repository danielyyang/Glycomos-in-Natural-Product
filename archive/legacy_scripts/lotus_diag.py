"""Quick diagnostic of LOTUS frozen.csv.gz"""
import pandas as pd

df = pd.read_csv("data/230106_frozen.csv.gz", dtype=str, low_memory=False)
print(f"Total rows: {len(df)}")
print(f"Unique InChIKeys: {df['structure_inchikey'].nunique()}")
print(f"Unique organisms: {df['organism_name'].nunique()}")
print()
print("Top 10 organisms:")
print(df["organism_name"].value_counts().head(10))
