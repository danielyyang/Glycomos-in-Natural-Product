import sys, os
from rdkit import Chem

SUGAR_SMILES_LIB = {
    "Glc": "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O", 
    "Gal": "OC[C@H]1O[C@@H](O)[C@@H](O)[C@@H](O)[C@@H]1O",
    "Man": "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O",
    "All": "OC[C@H]1O[C@@H](O)[C@H](O)[C@H](O)[C@@H]1O",   
}

def generate_wildcard_query(smiles):
    # Using SMARTS string manipulation
    sma = smiles
    
    # 1. Terminal O: ]O -> ]-[O,N,S,P]
    if sma.endswith("]O"):
        sma = sma[:-2] + "]-[O,N,S,P]"
        
    # 2. Exocyclic OH: (O) -> ([O,N,S,P])
    sma = sma.replace("(O)", "([O,N,S,P])")
    
    # 3. CH2OH group: OC[ -> [O,N,S,P]-[CH2]-[
    sma = sma.replace("OC[", "[O,N,S,P]-[CH2]-[")
    
    return Chem.MolFromSmarts(sma)

print("Testing SMARTS Wildcards")

queries = {}
for name, smi in SUGAR_SMILES_LIB.items():
    q = generate_wildcard_query(smi)
    if q: queries[name] = q
    else: print(f"Failed to generate query for {name}")

maltose_smi = "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O[C@H]1OC(CO)[C@@H](O)[C@H](O)[C@H]1O"
mol = Chem.MolFromSmiles(maltose_smi)
Chem.AssignStereochemistry(mol, force=True)

for name, q in queries.items():
    matches = mol.GetSubstructMatches(q, useChirality=True)
    if matches:
        print(f"Maltose matched {name}! ({len(matches)} times) At: {matches}")

cellobiose_smi = "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1OC1OC(CO)[C@@H](O)[C@H](O)[C@H]1O"
mol2 = Chem.MolFromSmiles(cellobiose_smi)
print("\nTesting Cellobiose")
for name, q in queries.items():
    matches = mol2.GetSubstructMatches(q, useChirality=True)
    if matches:
        print(f"Cellobiose matched {name}! ({len(matches)} times) At: {matches}")
