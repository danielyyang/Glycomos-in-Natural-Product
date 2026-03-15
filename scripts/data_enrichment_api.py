import pandas as pd
import requests
import time
import sys

sys.stdout.reconfigure(encoding='utf-8')

def fetch_isomeric_smiles_from_pubchem(inchi_key):
    """使用 InChIKey 向 PubChem 请求精确的 3D Isomeric SMILES"""
    if pd.isna(inchi_key) or not inchi_key: return None
    clean_key = str(inchi_key).replace("InChIKey=", "").strip()
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{clean_key}/property/IsomericSMILES/JSON"
    try:
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            return response.json()['PropertyTable']['Properties'][0]['IsomericSMILES']
    except Exception: pass
    return None

def format_dois(doi_string):
    if pd.isna(doi_string) or str(doi_string).strip() == "": return "No DOI found"
    return ", ".join([f"https://doi.org/{doi.strip()}" for doi in str(doi_string).split("|") if doi.strip()])

def run_enrichment_test(csv_path, limit=3):
    print("📡 启动 Phase 4 API 溯源与纠错雷达...")
    df = pd.read_csv(csv_path)

    # 重点抽查序列中含有 Hex 的行
    sample_df = df[df['Sugar_Sequence'].str.contains('Hex', na=False)].head(limit) if 'Sugar_Sequence' in df.columns else df.head(limit)
            
    for idx, row in sample_df.iterrows():
        print(f"\n" + "="*50)
        print(f"🆔 化合物 ID: {row['identifier']}")
        print(f"📝 常用名称: {row.get('name', 'N/A')}")
        print(f"🧬 原残缺 SMILES: {row['canonical_smiles']}")
        
        print("🔄 正在向 PubChem 呼叫高清三维结构...")
        iso_smiles = fetch_isomeric_smiles_from_pubchem(row['standard_inchi_key'])
        if iso_smiles and iso_smiles != row['canonical_smiles']:
            print(f"   ✅ 抢救成功! 精确结构: {iso_smiles}")
        else:
            print("   ❌ PubChem 无更优结构或网络超时。")
            
        print(f"📚 参考文献: {format_dois(row.get('dois', ''))}")
        print(f"🌿 生物来源 (Organism): {row.get('organisms', '未知')}")
        time.sleep(0.5) # 防止封禁

if __name__ == "__main__":
    # 指向最新的 CSV 输出文件
    run_enrichment_test("reports/Coconut_Sugar_Phase6.csv", limit=3)
