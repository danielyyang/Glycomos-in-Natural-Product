import sys, os
from rdkit import Chem
sys.path.append(os.path.abspath(os.path.dirname(__file__)))
from lib import sugar_sequence

sucrose = "OC[C@H]1O[C@@H](O[C@@]2(CO)O[C@H](CO)[C@@H](O)[C@@H]2O)[C@H](O)[C@@H](O)[C@H]1O"
mol = Chem.MolFromSmiles(sucrose)
ri = mol.GetRingInfo()
rings = [list(r) for r in ri.AtomRings()]

for i, ring in enumerate(rings):
    clean_mol, mapping = sugar_sequence.isolate_sugar_ring(mol, ring)
    mapped_ring = [mapping[idx] for idx in ring]
    sig = sugar_sequence.get_rs_signature_core(clean_mol, mapped_ring)
    
    print(f"Sucrose Ring {i} Size {len(ring)}:")
    print(f"  Isolated SMILES: {Chem.MolToSmiles(clean_mol)}")
    print(f"  Signature computed: {sig}")
    
    name = sugar_sequence.identify_monosaccharide_v2(mol, ring)
    print(f"  Identified Name: {name}")

