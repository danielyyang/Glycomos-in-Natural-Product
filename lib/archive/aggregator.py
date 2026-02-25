import requests
import pandas as pd
import logging
from typing import List, Dict, Optional
import time
import pubchempy as pcp
from urllib.error import HTTPError
import os

# 配置日志 (Configure logging)
logger = logging.getLogger(__name__)

class DataAggregator:
    """
    负责从多个数据源抓取糖脂数据 (Responsible for scraping glycolipid data from multiple sources)
    支持: GlyCosmos, COCONUT (API/Local), PubChem, KNApsack (Planned)
    """
    
    # GlyCosmos SPARQL Endpoint
    GLYCOSMOS_SPARQL_URL = "https://ts.glycosmos.org/sparql"
    # COCONUT API (Updated to v2 or v1 if available, fallback to simple search)
    COCONUT_API_URL = "https://coconut.naturalproducts.net/api/v2/search/simple" 

    def __init__(self):
        pass

    def fetch_glycosmos_data(self, output_file: str = "data/raw/glycosmos.csv") -> bool:
        """
        从 GlyCosmos 获取糖脂数据 (Fetch glycolipid data from GlyCosmos)
        """
        logging.info("Starting GlyCosmos SPARQL query...")
        
        query = """
        PREFIX glycan: <http://purl.jp/bio/12/glyco/glycan#>
        PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
        
        SELECT DISTINCT ?accession ?name ?smiles
        WHERE {
            ?structure a glycan:Glycan ;
                       glycan:has_accession ?accession .
            OPTIONAL { ?structure rdfs:label ?name } .
            OPTIONAL { ?structure glycan:has_smiles ?smiles } .
        }
        LIMIT 100
        """
        
        params = {
            "query": query,
            "format": "json"
        }
        
        try:
            response = requests.get(self.GLYCOSMOS_SPARQL_URL, params=params)
            response.raise_for_status()
            data = response.json()
            
            results = []
            bindings = data.get('results', {}).get('bindings', [])
            
            for item in bindings:
                entry = {
                    'Entry_ID': item.get('accession', {}).get('value'),
                    'Common_Name': item.get('name', {}).get('value', ''),
                    'SMILES': item.get('smiles', {}).get('value', ''),
                    'Source_DB': 'GlyCosmos'
                }
                if entry['SMILES']: 
                    results.append(entry)
            
            if results:
                df = pd.DataFrame(results)
                df.to_csv(output_file, index=False)
                logging.info(f"Successfully saved {len(df)} records to {output_file}")
                return True
            else:
                logging.warning("No results found from GlyCosmos.")
                return False
                
        except Exception as e:
            logger.error(f"Error scraping GlyCosmos: {e}")
            return False

    def fetch_coconut_data(self, keyword: str = "glycoside", limit: int = 50) -> List[Dict]:
        """
        从 COCONUT 提取数据 (Exctract data from COCONUT)
        """
        logging.info(f"Searching COCONUT for keyword: {keyword}...")
        results = []
        
        # Try finding the correct endpoint via simple string search
        search_url = f"{self.COCONUT_API_URL}?query={keyword}"
        
        try:
            # Note: The COCONUT API might be unstable or versioned.
            # If 404, we catch it.
            response = requests.get(search_url, timeout=15)
            response.raise_for_status()
            data = response.json()
            
            # Adjust parsing based on actual response structure
            items = data if isinstance(data, list) else data.get('results', [])
            if not items and 'compounds' in data:
                 items = data['compounds']

            for item in items[:limit]:
                # Adapt fields based on API response
                entry = {
                    'Entry_ID': item.get('coconut_id') or item.get('identifier', 'unknown'),
                    'Common_Name': item.get('name') or item.get('iupac_name', ''),
                    'SMILES': item.get('smiles') or item.get('canonical_smiles', ''),
                    'Source_DB': 'COCONUT',
                    'Category': keyword
                }
                if entry['SMILES']:
                    results.append(entry)
            
            logging.info(f"Fetched {len(results)} records from COCONUT.")
            
        except requests.exceptions.HTTPError as e:
            logger.error(f"COCONUT API Error ({e.response.status_code}): {e}")
            if e.response.status_code == 404:
                logger.warning("COCONUT API endpoint might have changed. Please verify documentation.")
        except Exception as e:
            logger.error(f"Error fetching from COCONUT: {e}")
            
        return results

    def fetch_pubchem_data(self, keyword: str = "glycoside", limit: int = 50) -> List[Dict]:
        """
        从 PubChem 获取数据 (Fetch data from PubChem)
        User recommended for 'Ginsenoside'.
        """
        logging.info(f"Searching PubChem for keyword: {keyword}...")
        results = []
        
        try:
            compounds = pcp.get_compounds(keyword, 'name', listkey_count=limit)
            
            for comp in compounds:
                entry = {
                    'Entry_ID': f"CID_{comp.cid}",
                    'Common_Name': keyword, # Ideally specific name, but keyword is safe
                    'SMILES': comp.isomeric_smiles,
                    'Source_DB': 'PubChem',
                    'Category': keyword
                }
                if entry['SMILES']:
                    results.append(entry)
                    
            logging.info(f"Fetched {len(results)} records from PubChem.")
            
        except Exception as e:
            logger.error(f"Error fetching from PubChem: {e}")
            
        return results

    def fetch_knapsack_data(self, keyword: str = "glycoside") -> List[Dict]:
        """
        KNApsack 数据获取占位符 (Placeholder for KNApsack fetcher)
        Currently implemented as a stub.
        """
        logging.info(f"Fetching from KNApsack (Stub) for {keyword}...")
        # TODO: Implement scraping for http://www.knapsackfamily.com/
        return []

    def fetch_from_local_file(self, file_path: str, source_name: str = "LocalFile", limit: int = 100) -> List[Dict]:
        """
        从本地文件加载数据 (Load data from local CSV/Excel)
        Suitable for downloaded databases like COCONUT.csv
        Expected columns: 'SMILES', (optional) 'Name'/'Common_Name', (optional) 'ID'/'Entry_ID'
        """
        logging.info(f"Loading local data from {file_path}...")
        results = []
        
        if not os.path.exists(file_path):
            logger.error(f"Local file not found: {file_path}")
            return results
            
        try:
            if file_path.endswith('.csv'):
                try:
                    df = pd.read_csv(file_path, encoding='utf-8-sig')
                except UnicodeDecodeError:
                    logger.warning("UTF-8 decode failed, trying GBK...")
                    df = pd.read_csv(file_path, encoding='gbk')
            elif file_path.endswith(('.xls', '.xlsx')):
                df = pd.read_excel(file_path)
            else:
                logger.error("Unsupported file format. Use CSV or Excel.")
                return results
                
            # Normalize columns
            cols = {c.lower(): c for c in df.columns}
            smiles_col = cols.get('smiles') or cols.get('canonical_smiles')
            if not smiles_col:
                logger.error("No 'SMILES' column found in local file.")
                return results
                
            name_col = cols.get('name') or cols.get('common_name') or cols.get('iupac_name')
            id_col = cols.get('id') or cols.get('entry_id') or cols.get('coconut_id')
            
            # Limit loading
            df = df.head(limit)
            
            for idx, row in df.iterrows():
                smiles = row[smiles_col]
                if pd.isna(smiles): continue
                
                entry = {
                    'Entry_ID': str(row[id_col]) if id_col else f"{source_name}_{idx}",
                    'Common_Name': str(row[name_col]) if name_col and not pd.isna(row[name_col]) else f"Unknown_{idx}",
                    'SMILES': smiles,
                    'Source_DB': source_name,
                    'Category': 'Local Import'
                }
                results.append(entry)
                
            logging.info(f"Loaded {len(results)} records from local file.")
            
        except Exception as e:
            logger.error(f"Error loading local file: {e}")
            
        return results

    def unify_results(self, data_list: List[List[Dict]], output_excel: str = "data/processed/unified_data.xlsx") -> pd.DataFrame:
        """
        数据清洗、去重并导出为 Excel (Data cleaning, deduplication, and export to Excel)
        """
        logging.info("Unifying and deduplicating results...")
        
        # Handle both flat list and list of lists
        if data_list and isinstance(data_list[0], list):
            all_entries = [item for sublist in data_list for item in sublist]
        else:
            all_entries = data_list
        
        if not all_entries:
            logging.warning("No data to unify.")
            return pd.DataFrame()
            
        df = pd.DataFrame(all_entries)
        
        # 1. 基础去重
        initial_count = len(df)
        if 'SMILES' in df.columns:
            df.drop_duplicates(subset=['SMILES'], keep='first', inplace=True)
            # 移除无效 SMILES
            df = df[df['SMILES'].str.len() > 5]
        
        logging.info(f"Deduplication complete: {initial_count} -> {len(df)} records.")
        
        # 3. 导出到 Excel
        try:
            os.makedirs(os.path.dirname(output_excel), exist_ok=True)
            df.to_excel(output_excel, index=False)
            logging.info(f"Unified data saved to {output_excel}")
        except Exception as e:
            logger.error(f"Error saving to Excel: {e}")
            
        return df

if __name__ == "__main__":
    aggregator = DataAggregator()
    # Test PubChem with Ginsenoside
    results = aggregator.fetch_pubchem_data("Ginsenoside", limit=5)
    print(results)
