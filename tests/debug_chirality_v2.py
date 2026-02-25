import sys
from rdkit import Chem

# Exact strings from generate_isomers.py
glc_s = "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O"
gal_s = "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O" 
man_s = "OC[C@H]1O[C@@H](O)[C@@H](O)[C@@H](O)[C@@H]1O"

REF_VARIANTS = {
    "Glc": [glc_s],
    "Gal": [gal_s],
    "Man": [man_s]
}

targets = {
    "Target_Glc": glc_s,
    "Target_Gal": gal_s,
    "Target_Man": man_s
}

for t_name, t_s in targets.items():
    t_mol = Chem.MolFromSmiles(t_s)
    print(f"\nChecking {t_name}...")
    for r_name, r_list in REF_VARIANTS.items():
        r_mol = Chem.MolFromSmiles(r_list[0])
        if t_mol.HasSubstructMatch(r_mol, useChirality=True):
            print(f"  MATCHES {r_name}")
        else:
            print(f"  NO MATCH {r_name}")
