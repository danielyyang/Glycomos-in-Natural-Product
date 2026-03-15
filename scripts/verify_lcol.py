import pandas as pd
from collections import Counter

df = pd.read_csv('reports/GlycoNP_Pipeline_Full_Cleaned.csv', low_memory=False)
seq = df['Sugar_Sequence'].dropna()
tokens = [t.strip() for s in seq for t in s.replace('\u2192', ',').replace('->', ',').split(',')]
c = Counter(tokens)

print('=== L-Col count ===')
print('L-Col:', c.get('L-Col', 0))
lcol_items = [(k, v) for k, v in c.items() if 'Col' in k]
print('All Col-related:', lcol_items[:20])
print()

print('=== Top 30 tokens ===')
for k, v in c.most_common(30):
    print(f'  {k:30s} {v:>8,}')
