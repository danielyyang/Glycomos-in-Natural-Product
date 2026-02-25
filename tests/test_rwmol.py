"""
Test RWMol logic to convert dummy atoms to hydrogens natively.
"""
import sys, os
from rdkit import Chem

cellobiose = "O[C@H]1[C@H](O)[C@@H](O)[C@H](O[C@@H]2[C@H](O)[C@@H](O)[C@H](O)[C@@H](CO)O2)[C@@H](CO)O1"
mol = Chem.MolFromSmiles(cellobiose)

def fix_dummy(frag_mol):
    rwmol = Chem.RWMol(frag_mol)
    for atom in rwmol.GetAtoms():
        if atom.GetAtomicNum() == 0:
            atom.SetAtomicNum(1) # H
            atom.SetIsotope(0)
    Chem.SanitizeMol(rwmol)
    return rwmol.GetMol()

# Fragment bond 18
bonds = [18] # Example
frag_mol = Chem.FragmentOnBonds(mol, bonds)
frags = Chem.GetMolFrags(frag_mol, asMols=True)

for i, frag in enumerate(frags):
    print(f"Frag {i}: {Chem.MolToSmiles(frag, isomericSmiles=True)}")
    clean = fix_dummy(frag)
    print(f"  Clean: {Chem.MolToSmiles(clean, isomericSmiles=True)}")
    Chem.AssignStereochemistry(clean, force=True, cleanIt=True)
    print(f"  Centers: {Chem.FindMolChiralCenters(clean)}")

