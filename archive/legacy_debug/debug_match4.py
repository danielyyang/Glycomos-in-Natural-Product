import sys, os
from rdkit import Chem
sys.path.append(os.path.abspath(os.path.dirname(__file__)))
from lib import sugar_sequence

smi = "OC[C@H]1O[C@@H](O[C@@]2(CO)O[C@H](CO)[C@@H](O)[C@@H]2O)[C@H](O)[C@@H](O)[C@H]1O"
mol = Chem.MolFromSmiles(smi)
ri = mol.GetRingInfo()
ring = [list(r) for r in ri.AtomRings()][1] # fructose
print(f"Target ring: {ring}")

clean_mol, mapping = sugar_sequence.isolate_sugar_ring(mol, ring)
matches1 = mol.GetSubstructMatches(clean_mol)
print(f"Total matches (default): {len(matches1)}")

matches2 = mol.GetSubstructMatches(clean_mol, useChirality=False)
print(f"Total matches (no chirality): {len(matches2)}")

# Could clean_mol have H explicitly?
has_h = any(a.GetSymbol() == 'H' for a in clean_mol.GetAtoms())
print(f"Has explicit H?: {has_h}")

# Could it be atomic number mismatch?
# What if clean_mol has a weird isotope?
has_iso = any(a.GetIsotope() != 0 for a in clean_mol.GetAtoms())
print(f"Has isotope?: {has_iso}")

# Let's see the SMILES
print(f"Sucrose SMILES: {Chem.MolToSmiles(mol)}")
print(f"Clean Mol SMILES: {Chem.MolToSmiles(clean_mol)}")

