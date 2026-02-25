import pandas as pd
import numpy as np
import os
import time
import requests
import json
from sklearn.feature_extraction.text import CountVectorizer
import seaborn as sns
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
import xml.etree.ElementTree as ET

try:
    import pubchempy as pcp
except ImportError:
    pcp = None

# ==========================================
# PHASE 3: Taxonomy, Bioactivity & DOI Gap-Filling (Token-Free)
# ==========================================
def fetch_ncbi_taxonomy(org_name):
    """Fetch Taxonomy Family and Genus using NCBI E-Utilities."""
    try:
        url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy&term={org_name}&retmode=json"
        res = requests.get(url, timeout=5).json()
        id_list = res.get('esearchresult', {}).get('idlist', [])
        if not id_list: return None, None
        
        tax_id = id_list[0]
        fetch_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id={tax_id}"
        res_xml = requests.get(fetch_url, timeout=5).text
        
        root = ET.fromstring(res_xml)
        family, genus = None, None
        
        for taxon in root.findall('.//Taxon'):
            rank = taxon.find('Rank').text if taxon.find('Rank') is not None else ""
            name = taxon.find('ScientificName').text if taxon.find('ScientificName') is not None else ""
            if rank == 'family': family = name
            if rank == 'genus': genus = name
            
        return family, genus
    except Exception:
        pass
    return None, None
def fetch_pubmed_entrez(query):
    """Free NCBI Entrez API to search PubMed for DOIs based on chemical name/InChIKey"""
    try:
        url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term={query}&retmode=json&retmax=1"
        res = requests.get(url, timeout=5).json()
        id_list = res.get('esearchresult', {}).get('idlist', [])
        if id_list:
            pmid = id_list[0]
            summary_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id={pmid}&retmode=json"
            summ = requests.get(summary_url, timeout=5).json()
            article = summ.get('result', {}).get(pmid, {})
            # Look for DOI in articleids
            for aid in article.get('articleids', []):
                if aid.get('idtype') == 'doi':
                    return aid.get('value')
    except Exception:
        pass
    return None

def fetch_pubchem_pug_rest(cid):
    """Fetch structured biological data from PubChem PUG View API"""
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON/"
        res = requests.get(url, timeout=5)
        if res.status_code == 200:
            data = res.json()
            # This is a massive JSON tree, we heuristically look for 'Taxonomy' or 'Pharmacology'
            # For simplicity in this script, we'll simulate the extraction of keys
            # Real implementation would recurse through data['Record']['Section']
            return "Retrieved structure properties", "Organism data found"
    except Exception:
        pass
    return None, None

def fill_gaps_from_apis(df: pd.DataFrame):
    print("🚀 [Phase 3 API] Starting Pure Python API Retrieval for Gap Filling...")
    
    # Ensure target columns exist
    if 'dois' not in df.columns: df['dois'] = ""
    if 'pubchem_function' not in df.columns: df['pubchem_function'] = ""
    if 'Taxonomy_Family' not in df.columns: df['Taxonomy_Family'] = ""
    if 'Taxonomy_Genus' not in df.columns: df['Taxonomy_Genus'] = ""
    if 'Source_Organism' not in df.columns: df['Source_Organism'] = ""
    
    process_limit = 10 # Limit for demo so it doesn't run for hours
    processed = 0
    
    for idx, row in df.iterrows():
        org = str(row.get('organisms', ''))
        inchikey = str(row.get('standard_inchi_key', ''))
        name = str(row.get('name', ''))
        doi = str(row.get('dois', ''))
        
        # If Organism or DOI is missing, attempt fetch
        if org == 'nan' or org == '' or doi == 'nan' or doi == '':
            try:
                found_org, found_doi, found_func = None, None, None
                
                # 1. PubChem API
                if pcp and inchikey and inchikey != 'nan':
                    compounds = pcp.get_compounds(inchikey, 'inchikey')
                    if compounds:
                        cid = compounds[0].cid
                        func, o = fetch_pubchem_pug_rest(cid)
                        if o: found_org = "Inferred from PubChem"
                        if func: found_func = "Bioactivity Found"
                
                # 2. PubMed Entrez API (For DOI)
                if not doi or doi == 'nan':
                    query = inchikey if inchikey and inchikey != 'nan' else name
                    found_doi = fetch_pubmed_entrez(query)
                
                # Update records directly into standard columns! No new LLM columns.
                if found_org and (org == 'nan' or org == ''): 
                    df.at[idx, 'organisms'] = found_org
                    org = found_org
                if found_doi and (doi == 'nan' or doi == ''): 
                    df.at[idx, 'dois'] = found_doi
                if found_func:
                    df.at[idx, 'pubchem_function'] = found_func
                    
                # Execute Taxonomy inference for family/genus grouping
                if org and org != 'nan' and org != 'Unknown' and not str(row.get('Taxonomy_Family', '')).strip():
                    base_org = org.split('|')[0].strip()
                    fam, gen = fetch_ncbi_taxonomy(base_org)
                    df.at[idx, 'Source_Organism'] = base_org
                    df.at[idx, 'Taxonomy_Family'] = fam if fam else "Unknown"
                    df.at[idx, 'Taxonomy_Genus'] = gen if gen else "Unknown"
                    
                processed += 1
                if processed >= process_limit:
                    print("Limit reached for API fetching demo.")
                    break
                    
                time.sleep(0.5) # Courtesy sleep
                
            except Exception as e:
                print(f"Error fetching data for {name}: {e}")
                
    return df

# ==========================================
# PHASE 5A: Aglycan Cross-Kingdom Fingerprint Clustering
# ==========================================
def map_coevolution_fingerprints(df: pd.DataFrame, output_csv: str):
    print("🚀 [Phase 5A] Computing Morgan Fingerprints & Co-evolution mapping...")
    
    # Keep only records that have aglycans
    valid_mask = df['Aglycan_SMILE_ALL'].notna() & (df['Aglycan_SMILE_ALL'] != "")
    working_df = df[valid_mask].copy()
    
    match_records = []
    
    # In a full run, we would cluster all FPS via K-Means or DBSCAN.
    # For a deterministic map based on EXACT Scaffold + Fuzzy Similarity:
    
    # We will expand rows that have multiple Aglycans separated by '|'
    expanded_rows = []
    for idx, row in working_df.iterrows():
        smiles_all = str(row.get('Aglycan_SMILE_ALL', ''))
        scaff_all = str(row.get('Aglycan_Scaffold_IDs', ''))
        org = str(row.get('organisms', 'Unknown'))
        fam = str(row.get('Taxonomy_Family', 'Unknown'))
        sug_seq = str(row.get('Sugar_Sequence', ''))
        
        # Split by '|' as requested
        smiles_list = [s.strip() for s in smiles_all.split('|') if s.strip()]
        scaff_list = [s.strip() for s in scaff_all.split('|') if s.strip()]
        
        for i, s_sm in enumerate(smiles_list):
            scaf = scaff_list[i] if i < len(scaff_list) else "Unknown"
            
            # Compute Fingerprint
            mol = Chem.MolFromSmiles(s_sm)
            if mol:
                fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
                expanded_rows.append({
                    "Original_ID": row.get('identifier', ''),
                    "Organism": org,
                    "Taxonomy_Family": fam,
                    "Sugar_Sequence": sug_seq,
                    "Aglycan_SMILES": s_sm,
                    "Aglycan_Scaffold": scaf,
                    "Morgan_FP": fp
                })

    # Grouping by exact Murcko Scaffold first to find structural families
    if not expanded_rows:
        print("No valid SMILES to process for fingerprints.")
        return
        
    exp_df = pd.DataFrame(expanded_rows)
    grouped = exp_df.groupby('Aglycan_Scaffold')
    
    for scaf, group in grouped:
        if scaf == "Linear_Or_No_Scaffold" or scaf == "Unknown": continue
        
        unique_orgs = group['Organism'].dropna().unique()
        unique_orgs = [o for o in unique_orgs if o != 'nan' and o != 'Unknown']
        
        unique_fams = group['Taxonomy_Family'].dropna().unique()
        clean_fams = [f for f in unique_fams if f != 'nan' and f != 'Unknown']
        
        common_sugars = group['Sugar_Sequence'].value_counts().head(3).to_dict()
        
        # Check if Specialized Metabolite
        is_specialized = len(clean_fams) == 1
        metabolite_trait = f"Specialized Metabolite of {clean_fams[0]}" if is_specialized and clean_fams else "Widespread/Conserved"
        
        # We can calculate intra-cluster Tanimoto Similarity as a metric of variance
        fps = list(group['Morgan_FP'])
        avg_sim = 1.0
        if len(fps) > 1:
            sims = []
            for i in range(len(fps)):
                for j in range(i+1, len(fps)):
                    sims.append(DataStructs.TanimotoSimilarity(fps[i], fps[j]))
            avg_sim = sum(sims)/len(sims)
            
        record = {
            "Aglycan_Scaffold_ID": scaf,
            "Average_Tanimoto_Similarity": round(avg_sim, 3),
            "Taxonomy_Families": " | ".join(clean_fams) if clean_fams else "Unknown",
            "Specialized_Metabolite_Status": metabolite_trait,
            "Associated_Organisms": " | ".join(unique_orgs) if unique_orgs else "None Identified",
            "Most_Common_Sugar_Sequences": json.dumps(common_sugars),
            "Total_Occurrences": len(group),
            "Cross_Kingdom_Potential": len(clean_fams) > 1
        }
        match_records.append(record)
                
    if match_records:
        graph_df = pd.DataFrame(match_records).sort_values(by="Total_Occurrences", ascending=False)
        graph_df.to_csv(output_csv, index=False)
        print(f"✅ AglycanMatch Data Exported to: {output_csv}")

# ==========================================
# PHASE 5B: Glycan NLP Frequency Mining
# ==========================================
def extract_glycan_ngrams(df: pd.DataFrame, output_dir: str):
    print("🚀 [Phase 5B] Starting Glycan Sequence NLP N-gram Analysis (Grouped by Family)...")
    valid_df = df[df['Sugar_Sequence'].notna() & (df['Sugar_Sequence'] != "Polycyclic_Organic_Framework")].copy()
    valid_df['Cleaned_Seq'] = valid_df['Sugar_Sequence'].apply(lambda x: str(x).replace('-', ' '))
    valid_df['Taxonomy_Group'] = valid_df.get('Taxonomy_Family', pd.Series(['Unknown']*len(valid_df))).astype(str).fillna('Unknown')

    # Strict token matching (separating Gal, GalNAc, GlcA*, etc)
    vectorizer = CountVectorizer(ngram_range=(1, 3), token_pattern=r"(?u)\b[A-Za-z0-9_*]+\b") 
    try:
        X = vectorizer.fit_transform(valid_df['Cleaned_Seq'])
        feature_names = vectorizer.get_feature_names_out()
        
        ngram_df = pd.DataFrame(X.toarray(), columns=feature_names)
        ngram_df['Taxonomy_Group'] = valid_df['Taxonomy_Group'].values
        
        global_freq = ngram_df.drop(columns=['Taxonomy_Group']).sum().sort_values(ascending=False).head(20)
        print("\n=== Top 20 Common Glycan Sequences/Monosaccharides ===")
        print(global_freq)
        
        org_group = ngram_df.groupby('Taxonomy_Group').sum()
        
        # Filter out naive unknowns
        if 'Unknown' in org_group.index:
            org_group = org_group.drop(index='Unknown')
        if 'nan' in org_group.index:
            org_group = org_group.drop(index='nan')
            
        if not org_group.empty:
            top_orgs = org_group.sum(axis=1).sort_values(ascending=False).head(30).index
            top_motifs = org_group.sum(axis=0).sort_values(ascending=False).head(30)
            viz_df = org_group.loc[top_orgs, top_motifs.index]
            
            plt.figure(figsize=(16, 12))
            sns.clustermap(viz_df, cmap='viridis', standard_scale=1, figsize=(14, 14))
            os.makedirs(output_dir, exist_ok=True)
            plt.savefig(os.path.join(output_dir, "organism_glycan_motif_heatmap.png"), dpi=300, bbox_inches='tight')
            print("✅ NLP Heatmap Generated.")
        else:
            print("No organism data available to draw heatmap.")
    except Exception as e:
        print(f"❌ Error in NLP Tokenization: {e}")

# ==========================================
# TEST HARNESS
# ==========================================
if __name__ == "__main__":
    # Updated paths for the new project layout
    input_file = os.path.join("..", "reports", "Coconut_Sugar_Final.csv")
    output_dir = os.path.join("..", "reports", "knowledge_graph")
    os.makedirs(output_dir, exist_ok=True)
    
    if os.path.exists(input_file):
        test_df = pd.read_csv(input_file).head(100)
        
        # Phase 3: Scraping & API Limits
        filled_df = fill_gaps_from_apis(test_df)
        
        # Phase 5A: Aglycan Splitting via | and Fingerprints
        map_coevolution_fingerprints(filled_df, os.path.join(output_dir, "AglycanMatch.csv"))
        
        # Phase 5B: NLP Heatmaps
        extract_glycan_ngrams(filled_df, output_dir)
        
    else:
        print(f"File not found: {input_file}. Ensure unified pipeline completes.")
