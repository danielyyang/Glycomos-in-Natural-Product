import sys, os
from rdkit import Chem
sys.path.append(os.path.abspath(os.path.dirname(__file__)))
from lib import sugar_sequence

maltose = "O[C@H]1[C@H](O)[C@@H](O)[C@H](O[C@H]2[C@H](O)[C@@H](O)[C@H](O)[C@@H](CO)O2)[C@@H](CO)O1"
mol = Chem.MolFromSmiles(maltose)
ri = mol.GetRingInfo()
rings = [list(r) for r in ri.AtomRings()]

for i, ring in enumerate(rings):
    clean_mol, mapping = sugar_sequence.isolate_sugar_ring(mol, ring)
    mapped_ring = [mapping[idx] for idx in ring]
    sig = sugar_sequence.get_rs_signature_core(clean_mol, mapped_ring)
    rs_key = (len(mapped_ring), sig)
    print(f"Ring {i}: Size {len(mapped_ring)}, Signature: {sig}")
    print(f"  -> SMILES: {Chem.MolToSmiles(clean_mol)}")
    print(f"  -> Match in Lib? {rs_key in sugar_sequence.RS_LIBRARY}")
    if rs_key in sugar_sequence.RS_LIBRARY:
        print(f" -> {sugar_sequence.RS_LIBRARY[rs_key]}")

