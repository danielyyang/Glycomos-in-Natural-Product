import sys, os
from rdkit import Chem
sys.path.append(os.path.abspath(os.path.dirname(__file__)))
from lib import sugar_sequence

maltose = "O[C@H]1[C@H](O)[C@@H](O)[C@H](O[C@H]2[C@H](O)[C@@H](O)[C@H](O)[C@@H](CO)O2)[C@@H](CO)O1"
mol = Chem.MolFromSmiles(maltose)
ri = mol.GetRingInfo()
rings = [list(r) for r in ri.AtomRings()]
ring = rings[0] # Acceptor

clean_mol, mapping = sugar_sequence.isolate_sugar_ring(mol, ring)
# clean_mol already has stereocenters assigned
centers = Chem.FindMolChiralCenters(clean_mol, includeUnassigned=True)
print(f"Centers (clean_mol): {centers}")

mapped_ring = [mapping[idx] for idx in ring]
print(f"Mapped ring: {mapped_ring}")
for idx in mapped_ring:
    if idx != -1: print(f"  Atom {idx}: {clean_mol.GetAtomWithIdx(idx).GetSymbol()}")

sig = sugar_sequence.get_rs_signature_core(clean_mol, mapped_ring)
print(f"Signature: {sig}")

from lib import sugar_utils
u = sugar_utils.find_mapped_sugar_units(mol)[0]
print(f"Algorithm Unit 0 name: {u['name']}")

