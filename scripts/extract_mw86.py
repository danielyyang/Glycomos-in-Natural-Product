"""
Extract 10 Unknown_Sugar_MW_86 SMILES for investigation.
提取 10 个 Unknown_Sugar_MW_86 的 SMILES 用于排查切割伪影。
"""
import pandas as pd

df = pd.read_csv('reports/GlycoNP_Pipeline_Full_Cleaned.csv', low_memory=False)
mw86 = df[df['Sugar_Sequence'].str.contains('Unknown_Sugar_MW_86', na=False)]
print(f'Total rows with MW_86: {len(mw86)}')
print()
print('=== 10 Unknown_Sugar_MW_86 Sample SMILES ===')
print()
for i, (_, row) in enumerate(mw86.head(10).iterrows()):
    cid = str(row.get('identifier', '?'))[:20]
    glySmiles = str(row.get('Glycan_SMILES', ''))[:120]
    fullSmiles = str(row.get('canonical_smiles', ''))[:120]
    seq = str(row.get('Sugar_Sequence', ''))[:80]
    print(f'{i+1}. ID={cid}')
    print(f'   Seq={seq}')
    print(f'   Glycan={glySmiles}')
    print(f'   Full={fullSmiles}')
    print()
