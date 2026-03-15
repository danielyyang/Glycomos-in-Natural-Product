import sys, os
from rdkit import Chem
sys.path.append(os.path.abspath(os.path.dirname(__file__)))
from lib import sugar_sequence

maltose = "O[C@H]1[C@H](O)[C@@H](O)[C@H](O[C@H]2[C@H](O)[C@@H](O)[C@H](O)[C@@H](CO)O2)[C@@H](CO)O1"
mol = Chem.MolFromSmiles(maltose)
ri = mol.GetRingInfo()
rings = [list(r) for r in ri.AtomRings()]

for i, ring in enumerate(rings):
    name = sugar_sequence.identify_monosaccharide_v2(mol, ring)
    print(f"Ring {i} -> {name}")

