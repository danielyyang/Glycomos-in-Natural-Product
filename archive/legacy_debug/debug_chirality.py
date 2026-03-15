import sys
import os
from rdkit import Chem
from rdkit.Chem import AssignStereochemistry

# Refs from process_sugar.py
REF_VARIANTS = {
    "Glc": [
        "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O", # beta
        "OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O",  # alpha
    ],
    "Gal": [
        "OC[C@H]1O[C@@H](O)[C@@H](O)[C@@H](O)[C@@H]1O", # beta
        "OC[C@H]1O[C@H](O)[C@@H](O)[C@@H](O)[C@@H]1O",  # alpha
    ],
     "Man": [
        "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O", # beta-D-Man
        "OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@H]1O",  # alpha-D-Man
    ]
}

def check(target_smiles, name):
    target = Chem.MolFromSmiles(target_smiles)
    if not target:
        print(f"Failed to load target {name}")
        return
    
    print(f"\nChecking {name} ({target_smiles})")
    
    for ref_name, variants in REF_VARIANTS.items():
        for i, v_smiles in enumerate(variants):
            ref = Chem.MolFromSmiles(v_smiles)
            if target.HasSubstructMatch(ref, useChirality=True):
                print(f"  MATCHES {ref_name} (Variant {i})")
            else:
                pass
                # print(f"  NO MATCH {ref_name} (Variant {i})")

if __name__ == "__main__":
    # Test cases:
    glc_beta = "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O"
    check(glc_beta, "Beta-D-Glc")
    
    gal_beta = "OC[C@H]1O[C@@H](O)[C@@H](O)[C@@H](O)[C@@H]1O"
    check(gal_beta, "Beta-D-Gal")
    
    man_beta = "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O"
    check(man_beta, "Beta-D-Man")
