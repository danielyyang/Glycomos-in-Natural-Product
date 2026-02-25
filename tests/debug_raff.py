import sys, os
from rdkit import Chem
sys.path.append(os.path.abspath(os.path.dirname(os.path.dirname(__file__))))
from lib import sugar_utils, sugar_sequence

raffinose = "C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@@]2([C@H]([C@@H]([C@H](O2)CO)O)O)CO)O)O)O)O[C@H]3[C@@H]([C@H]([C@@H]([C@H](O3)CO)O)O)O"
mol = Chem.MolFromSmiles(raffinose)

for ri, ring in enumerate(mol.GetRingInfo().AtomRings()):
     clean_mol, mapping = sugar_sequence.isolate_sugar_ring(mol, list(ring))
     if clean_mol is None:
         print(f"Ring {ri}: Failed to isolate")
         continue
     mapped_ring = [mapping[idx] for idx in ring]
     sig = sugar_sequence.get_rs_signature_core(clean_mol, mapped_ring)
     name = sugar_sequence.identify_monosaccharide_v2(mol, ring)
     print(f"Ring {ri}: size {len(ring)} -> sig {sig} -> name {name}")
