import pandas as pd
try:
    df = pd.read_csv('data/processed/Coconut_Sugar_Sort_Enriched_v2.csv')
    pd.set_option('display.max_columns', None)
    pd.set_option('display.max_colwidth', 100)
    print(df[['Glycan_SMILES', 'Refined_Sequence']].head(20).to_string())
except Exception as e:
    print(e)
