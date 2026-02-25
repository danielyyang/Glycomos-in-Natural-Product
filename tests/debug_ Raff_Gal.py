import sys, os
from rdkit import Chem
sys.path.append(os.path.abspath(os.path.dirname(os.path.dirname(__file__))))
from lib import sugar_sequence, sugar_utils

gal_smiles = "C([C@@H]1[C@@H]([C@@H]([C@H]([C@H](O1)O)O)O)O)O" # User's pubchem Gal
raff_smiles = "C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@@]2([C@H]([C@@H]([C@H](O2)CO)O)O)CO)O)O)O)O[C@H]3[C@@H]([C@H]([C@@H]([C@H](O3)CO)O)O)O"

mol_gal = Chem.MolFromSmiles(gal_smiles)
mol_raff = Chem.MolFromSmiles(raff_smiles)

print("Galactose PubChem:")
for ri, ring in enumerate(mol_gal.GetRingInfo().AtomRings()):
     clean_mol, mapping = sugar_sequence.isolate_sugar_ring(mol_gal, list(ring))
     mapped_ring = [mapping[idx] for idx in ring]
     sig = sugar_sequence.get_rs_signature_core(clean_mol, mapped_ring)
     name = sugar_sequence.identify_monosaccharide_v2(mol_gal, ring)
     print(f"  Ring {ri}: size {len(ring)} -> sig {sig} -> name {name}")

print("\nRaffinose Test Case:")
for ri, ring in enumerate(mol_raff.GetRingInfo().AtomRings()):
     clean_mol, mapping = sugar_sequence.isolate_sugar_ring(mol_raff, list(ring))
     mapped_ring = [mapping[idx] for idx in ring]
     sig = sugar_sequence.get_rs_signature_core(clean_mol, mapped_ring)
     name = sugar_sequence.identify_monosaccharide_v2(mol_raff, ring)
     print(f"  Ring {ri}: size {len(ring)} -> sig {sig} -> name {name}")
