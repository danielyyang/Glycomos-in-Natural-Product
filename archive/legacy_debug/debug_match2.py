import sys, os
from rdkit import Chem
sys.path.append(os.path.abspath(os.path.dirname(__file__)))
from lib import sugar_sequence

smi = "OC[C@H]1O[C@@H](O[C@@]2(CO)O[C@H](CO)[C@@H](O)[C@@H]2O)[C@H](O)[C@@H](O)[C@H]1O"
mol = Chem.MolFromSmiles(smi)
ri = mol.GetRingInfo()
ring = [list(r) for r in ri.AtomRings()][1] # fructose
print(f"Target ring: {ring}")
for idx in ring:
    print(f"  Atom {idx} in mol is {mol.GetAtomWithIdx(idx).GetSymbol()}")
    
clean_mol, mapping = sugar_sequence.isolate_sugar_ring(mol, ring)
match = mol.GetSubstructMatch(clean_mol)
clean_to_target = {c: t for c, t in enumerate(match)}
target_to_clean = {t: c for c, t in clean_to_target.items()}

print("\nMatch mapping:")
for c, t in clean_to_target.items():
    print(f"  clean[{c}] ({clean_mol.GetAtomWithIdx(c).GetSymbol()}) -> mol[{t}] ({mol.GetAtomWithIdx(t).GetSymbol()})")

mapped_ring = [target_to_clean.get(old_r, -1) for old_r in ring]
print(f"\nMapped ring: {mapped_ring}")
for idx in mapped_ring:
    print(f"  Atom {idx} in clean_mol is {clean_mol.GetAtomWithIdx(idx).GetSymbol()}")
