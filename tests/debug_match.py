import sys, os
from rdkit import Chem
sys.path.append(os.path.abspath(os.path.dirname(__file__)))
from lib import sugar_sequence

smi = "OC[C@H]1O[C@@H](O[C@@]2(CO)O[C@H](CO)[C@@H](O)[C@@H]2O)[C@H](O)[C@@H](O)[C@H]1O"
mol = Chem.MolFromSmiles(smi)
ri = mol.GetRingInfo()
ring = [list(r) for r in ri.AtomRings()][1] # fructose

clean_mol, mapping = sugar_sequence.isolate_sugar_ring(mol, ring)
match = mol.GetSubstructMatch(clean_mol)
print(f"Match found? {bool(match)}")

# If we use `useChirality=False` explicitly?
match_no_chi = mol.GetSubstructMatch(clean_mol, useChirality=False)
print(f"Match without chirality? {bool(match_no_chi)}")

# What if clean_mol is failing because of the new Hydrogens?
# clean_mol has explicit OH where it used to be a glycosidic bond. SubstructMatch is exact.
# In Sucrose, the glycosidic bond is between Glc C1 and Fru C2.
# So Fru C2 in Sucrose is connected to Glc O1.
# In clean_mol, Fru C2 is connected to an OH! 
# But `mol` (Sucrose) doesn't have an OH on Fru C2, it has an O connected to Glc C1!
# Wait! In RDKit, `GetSubstructMatch` matches atoms by element.
# clean_mol has an OH. mol has an O-C. 
# They DO NOT match perfectly if we use standard molecules.
# Oh! The dummy replacement made it an OH.
