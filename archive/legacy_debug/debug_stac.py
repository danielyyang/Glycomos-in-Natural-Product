import sys, os
from rdkit import Chem
sys.path.append(os.path.abspath(os.path.dirname(__file__)))
from lib import sugar_sequence

# Stachyose (Gal-a1,6-Gal-a1,6-Glc-a1,2-Fru)
stachyose = "OC[C@H]1O[C@H](O[C@H]2[C@H](O)[C@@H](O)[C@H](O[C@@]3(CO)O[C@H](CO)[C@@H](O)[C@@H]3O)[C@H](CO[C@H]4O[C@H](CO[C@H]5O[C@H](CO)[C@H](O)[C@H](O)[C@H]5O)[C@H](O)[C@H](O)[C@H]4O)O2)[C@H](O)[C@@H](O)[C@@H]1O"
mol = Chem.MolFromSmiles(stachyose)
ri = mol.GetRingInfo()
rings = [list(r) for r in ri.AtomRings()]

for i, ring in enumerate(rings):
    clean_mol, mapping = sugar_sequence.isolate_sugar_ring(mol, ring)
    mapped_ring = [mapping[idx] for idx in ring]
    sig = sugar_sequence.get_rs_signature_core(clean_mol, mapped_ring)
    rs_key = (len(mapped_ring), sig)
    name = sugar_sequence.identify_monosaccharide_v2(mol, ring)
    
    print(f"Stachyose Ring {i}:")
    print(f"  Size {len(mapped_ring)}, Signature: {sig}")
    print(f"  Identified Name: {name}")
    print(f"  RS_LIBRARY hit: {rs_key in sugar_sequence.RS_LIBRARY}")
    if rs_key in sugar_sequence.RS_LIBRARY:
         print(f"  -> {sugar_sequence.RS_LIBRARY[rs_key]}")

