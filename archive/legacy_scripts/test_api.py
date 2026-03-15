import requests
import urllib.parse
import json

def test_iupac_to_organism(iupac_name):
    print(f"Testing IUPAC: {iupac_name}")
    # 1. Get CID from IUPAC
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{urllib.parse.quote(iupac_name)}/cids/JSON"
    res = requests.get(url)
    if res.status_code != 200:
        print("Failed to get CID")
        return
        
    cids = res.json().get('IdentifierList', {}).get('CID', [])
    if not cids:
        print("No CIDs found")
        return
        
    cid = cids[0]
    print(f"Found CID: {cid}")
    
    # 2. Get Taxonomy data from PUG View
    tax_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON?heading=Taxonomy"
    tax_res = requests.get(tax_url)
    if tax_res.status_code != 200:
        print("Failed to get Taxonomy page. Might not have taxonomy info.")
        return
        
    data = tax_res.json()
    record = data.get('Record', {})
    
    # PubChem Taxonomy data is deeply nested. Let's try to extract organism names.
    organisms = set()
    try:
        sections = record.get('Section', [])
        for sec in sections:
            if sec.get('TOCHeading') == 'Taxonomy':
                for info in sec.get('Information', []):
                    # Usually Organism is in the Name or Value
                    if 'Name' in info and 'Value' in info:
                        val = info['Value']['StringWithMarkup'][0]['String']
                        organisms.add(val)
    except Exception as e:
        print("Error parsing taxonomy:", e)
        
    print("Organisms found in PubChem:", list(organisms))

def test_organism_to_family(organism):
    print(f"\nTesting Organism: {organism}")
    url = f"https://api.gbif.org/v1/species/match?name={urllib.parse.quote(organism)}"
    res = requests.get(url)
    if res.status_code == 200:
        data = res.json()
        print(f"Family: {data.get('family', 'Not Found')}")
    else:
        print("GBIF API failed")

if __name__ == "__main__":
    # Test with a known common IUPAC name or synonym that has taxonomy
    # E.g. "Ginsenoside Rb1" or "beta-D-Glucopyranose"
    test_iupac_to_organism("Ginsenoside Rb1")
    test_organism_to_family("Panax ginseng")
    test_organism_to_family("Oryza sativa")
