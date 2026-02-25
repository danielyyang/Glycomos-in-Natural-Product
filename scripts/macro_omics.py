import pandas as pd
import numpy as np
import os
import re
import argparse
import warnings
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from sklearn.cluster import DBSCAN
from sklearn.feature_extraction.text import CountVectorizer
import matplotlib.pyplot as plt
import seaborn as sns

warnings.filterwarnings('ignore')

def fetch_taxonomy(organism_name):
    """
    Simulated fast taxonomy resolution. 
    In production, this could query NCBI Taxonomy API.
    """
    if pd.isna(organism_name) or str(organism_name).strip() == "":
        return "Unknown_Family", "Unknown"
        
    genus = str(organism_name).split()[0].capitalize()
    return f"{genus}aceae", genus

def custom_sugar_tokenizer(text):
    """Splits 'Glc(1-4)GalNAc(1-2)Rha' -> ['Glc', 'GalNAc', 'Rha']"""
    if not isinstance(text, str):
        return []
    clean_text = re.sub(r'\([0-9\?\-a-zA-Z\alpha\beta\,]+\)', ' ', text)
    tokens = [t.strip() for t in clean_text.split() if t.strip()]
    return tokens

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--limit", type=int, default=0, help="Limit rows processing")
    args = parser.parse_args()

    base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    enriched_csv = os.path.join(base_dir, "reports", "Coconut_Sugar_Enriched.csv")
    final_csv = os.path.join(base_dir, "reports", "Coconut_Sugar_Final.csv")
    
    input_csv = enriched_csv if os.path.exists(enriched_csv) else final_csv
    if not os.path.exists(input_csv):
        print("Data files not found. Please run unified_pipeline.py first.")
        return
        
    df = pd.read_csv(input_csv, low_memory=False)
    
    if args.limit > 0:
        df = df.head(args.limit).copy()
        print(f"Limiting to {args.limit} rows.")
        
    print("\n--- Phase 5.0: Taxonomy Resolution ---")
    families = []
    genera = []
    for org in tqdm(df['organisms'], desc="Resolving Taxonomy"):
        family, genus = fetch_taxonomy(org)
        families.append(family)
        genera.append(genus)
    df['Taxonomy_Family'] = families
    df['Taxonomy_Genus'] = genera

    print("\n--- Phase 5: Per-Class Combinatorial Macro-omics Analysis ---")
    
    classes = df['chemical_super_class'].dropna().unique()
    
    for cls in classes:
        if cls == 'No Classified':
            continue
            
        print(f"\nProcessing Class >> {cls} <<")
        class_df = df[df['chemical_super_class'] == cls].copy()
        
        if class_df.empty:
            continue
            
        # Create output directory for this specific class
        safe_cls_name = str(cls).replace('/', '_').replace('\\', '_').replace(' ', '_').replace('>', '')
        cls_dir = os.path.join(base_dir, "reports", safe_cls_name)
        os.makedirs(cls_dir, exist_ok=True)
        
        # ----------------------------------------------------
        # 1. Aglycan Clustering within this class
        # ----------------------------------------------------
        aglycan_smiles_list = class_df['Aglycan_SMILE_ALL'].dropna().tolist()
        unique_aglycans = set()
        for ags in aglycan_smiles_list:
            for s in str(ags).split('|'):
                s = s.strip()
                if s and s != 'nan':
                    unique_aglycans.add(s)
                    
        unique_aglycans = list(unique_aglycans)
        cluster_mapping = {}
        
        if unique_aglycans:
            fps = []
            valid_smiles = []
            for s in unique_aglycans:
                try:
                    mol = Chem.MolFromSmiles(s)
                    if mol:
                        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
                        fps.append(fp)
                        valid_smiles.append(s)
                except:
                    pass
            
            n = len(fps)
            if n > 1:
                dist_matrix = np.zeros((n, n))
                for i in range(n):
                    for j in range(i, n):
                        sim = DataStructs.TanimotoSimilarity(fps[i], fps[j])
                        dist = 1.0 - sim
                        dist_matrix[i, j] = dist
                        dist_matrix[j, i] = dist
                
                # DBSCAN Cluster (Threshold = 0.8 Sim -> eps=0.2)
                db = DBSCAN(eps=0.2, min_samples=2, metric='precomputed')
                labels = db.fit_predict(dist_matrix)
                
                for s, l in zip(valid_smiles, labels):
                    cluster_mapping[s] = f"Cluster_{str(l).zfill(2)}" if l != -1 else "Orphan"
            elif n == 1:
                cluster_mapping[valid_smiles[0]] = "Cluster_00"
                
        def assign_cluster(ag_str):
            if pd.isna(ag_str): return "None"
            parts = [s.strip() for s in str(ag_str).split('|') if s.strip()]
            if not parts: return "None"
            return cluster_mapping.get(parts[0], "Unknown")
            
        class_df['Aglycan_Cluster'] = class_df['Aglycan_SMILE_ALL'].apply(assign_cluster)
        
        # Save cluster definitions for this class
        cluster_out = os.path.join(cls_dir, f"{safe_cls_name}_Aglycan_Clusters.csv")
        cluster_data = pd.DataFrame(list(cluster_mapping.items()), columns=['Aglycan_SMILES', 'Cluster_ID']).sort_values('Cluster_ID')
        cluster_data.to_csv(cluster_out, index=False)
        
        # ----------------------------------------------------
        # 2. NLP Sequence Decomposition & Heatmap Generation
        # ----------------------------------------------------
        valid_seqs = class_df['Sugar_Sequence'].dropna()
        valid_seqs = valid_seqs[valid_seqs != ""]
        
        if valid_seqs.empty:
            print(f"  [!] No valid sugar sequences found in {cls}. Skipping heatmaps.")
            continue
            
        vectorizer = CountVectorizer(ngram_range=(1, 3), tokenizer=custom_sugar_tokenizer, min_df=1)
        try:
            X = vectorizer.fit_transform(valid_seqs)
            vocab = vectorizer.get_feature_names_out()
        except ValueError:
            print(f"  [!] Could not vectorize sequences for {cls}. Skipping.")
            continue
            
        seq_df = class_df.loc[valid_seqs.index].copy()
        
        # Filter top N most frequent N-grams across this class to prevent huge matrices
        # Yes, it automatically selects high frequency motifs (Top 30)!
        top_n = min(30, len(vocab))
        col_sums = X.sum(axis=0).A1
        top_indices = col_sums.argsort()[::-1][:top_n]
        top_vocab = vocab[top_indices]
        
        # Dense matrix for top N motifs
        X_top = X[:, top_indices].toarray()
        motif_df = pd.DataFrame(X_top, columns=top_vocab, index=seq_df.index)
        
        def plot_crosstab_heatmap(group_col, title, filename):
            # Group by specified column and sum motif frequencies
            seq_df[group_col] = seq_df[group_col].fillna("Unknown")
            grouped = motif_df.groupby(seq_df[group_col]).sum()
            
            # Prune empty rows/cols
            grouped = grouped.loc[(grouped.sum(axis=1) > 0), (grouped.sum(axis=0) > 0)]
            
            if grouped.shape[0] < 2 or grouped.shape[1] < 2:
                print(f"  [~] Not enough data variance to plot {filename}. Requires at least 2x2 matrix.")
                return
                
            plt.figure(figsize=(10, 8))
            try:
                sns.clustermap(
                    grouped, 
                    cmap="YlGnBu", 
                    standard_scale=1, # Scale across motif columns to highlight preference
                    figsize=(14, 10),
                    linewidths=.5,
                    annot=False
                )
                plt.title(title, pad=20, fontsize=16)
                plt.savefig(os.path.join(cls_dir, filename), dpi=300, bbox_inches='tight')
                plt.close('all')
                print(f"  --> Saved {filename}")
            except Exception as e:
                print(f"  [Error] Failed to plot {filename}: {e}")

        # Generate the 3 Combinatorial Heatmaps
        print(f"  Generating NLP Clustered Heatmaps (Top {top_n} Motifs)...")
        plot_crosstab_heatmap('Aglycan_Cluster', f'{cls}\nAglycan Cluster vs Sugar Motifs', 'Heatmap_Aglycan_vs_Sugar.png')
        plot_crosstab_heatmap('organisms', f'{cls}\nOrganism vs Sugar Motifs', 'Heatmap_Organism_vs_Sugar.png')
        plot_crosstab_heatmap('Taxonomy_Family', f'{cls}\nTaxonomy Family vs Sugar Motifs', 'Heatmap_Family_vs_Sugar.png')
        
    print("\nAll macro-omics processing complete! Please check the reports/ directories.")

if __name__ == "__main__":
    main()
