import pandas as pd
df = pd.read_csv('data/processed/Coconut_Sugar_Sort_Enriched.csv')
print(df[['Refined_Sequence']].head(20).to_string())
