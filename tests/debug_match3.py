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
matches = mol.GetSubstructMatches(clean_mol)
print(f"Total matches found: {len(matches)}")

for i, match in enumerate(matches):
    overlap = sum(1 for t in match if t in ring)
    overlap_list = [t for t in match if t in ring]
    print(f"Match {i}: {match}")
    print(f"  Overlap count: {overlap} / {len(ring)}")
    print(f"  Overlap atoms: {overlap_list}")

