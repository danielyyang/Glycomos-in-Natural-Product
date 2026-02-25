import sys, os
from rdkit import Chem
sys.path.append(os.path.abspath(os.path.dirname(os.path.dirname(__file__))))
from lib import sugar_sequence

smiles = "O[C@H]1[C@H](O)[C@@H](O)[C@H](O[C@H]2[C@H](O)[C@@H](O)[C@H](O)[C@@H](CO)O2)[C@@H](CO)O1" # Maltose
mol = Chem.MolFromSmiles(smiles)

print("Maltose Rings:")
for ri, ring in enumerate(mol.GetRingInfo().AtomRings()):
     clean_mol, mapping = sugar_sequence.isolate_sugar_ring(mol, list(ring))
     canon_smi = Chem.MolToSmiles(clean_mol, isomericSmiles=True)
     flat_smi = Chem.MolToSmiles(clean_mol, isomericSmiles=False)
     print(f"Ring {ri}: Canonical: {canon_smi} | Flat: {flat_smi}")
