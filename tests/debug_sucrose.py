import traceback
from rdkit import Chem
import sys, os
sys.path.append(os.path.abspath(os.path.dirname(__file__)))
from lib import sugar_utils
from lib import sugar_sequence

print("Loading Sucrose...")
smi = "OC[C@H]1O[C@@H](O[C@@]2(CO)O[C@H](CO)[C@@H](O)[C@@H]2O)[C@H](O)[C@@H](O)[C@H]1O"
mol = Chem.MolFromSmiles(smi)

try:
    print("Running get_sugar_units...")
    units, atom_map = sugar_utils.get_sugar_units(mol)
    print("Success!")
except Exception as e:
    print(f"Crash in get_sugar_units: {e}")
    traceback.print_exc()

ri = mol.GetRingInfo()
valid_rings = []
for ring in ri.AtomRings():
    atoms = [mol.GetAtomWithIdx(i) for i in ring]
    if [a.GetSymbol() for a in atoms].count("O") == 1:
        valid_rings.append(list(ring))

from lib.sugar_sequence import get_rs_signature_core, isolate_sugar_ring, RS_LIBRARY

for sid, ring in enumerate(valid_rings):
    print(f"\n--- Ring {sid} ---")
    try:
        clean, mapping = isolate_sugar_ring(mol, list(ring))
        if clean and mapping:
            Chem.AssignStereochemistry(clean, force=True, cleanIt=True)
            mapped_ring = [mapping[idx] for idx in ring]
            print(f"Mapped ring: {mapped_ring}")
            
            # Now let's trace the signature generation step-by-step
            o_idx = None
            for idx in mapped_ring:
                if clean.GetAtomWithIdx(idx).GetSymbol() == "O": o_idx=idx; break
                
            o_atom = clean.GetAtomWithIdx(o_idx)
            neighbors = [n for n in o_atom.GetNeighbors() if n.GetIdx() in mapped_ring]
            if len(neighbors) != 2: 
                 print("Not 2 neighbors to O")
                 continue
            n1, n2 = neighbors
            
            def score_c1(atom_idx):
                atom = clean.GetAtomWithIdx(atom_idx)
                score = 0
                for n in atom.GetNeighbors():
                    if n.GetIdx() == o_idx: continue
                    sym = n.GetSymbol()
                    if sym in ["O", "N", "S", "P"]: score += 2
                    if sym == "C": score += 0.5 
                return score
                
            s1 = score_c1(n1.GetIdx())
            s2 = score_c1(n2.GetIdx())
            print(f"n1={n1.GetIdx()} score={s1}, n2={n2.GetIdx()} score={s2}")
            
            sig = get_rs_signature_core(clean, mapped_ring)
            print(f"Signature: {sig}")
            if sig in RS_LIBRARY:
                print(f"Matched: {RS_LIBRARY[sig]}")
            else:
                print("NO MATCH")
    except Exception as e:
        print(f"Crash isolating or scoring ring {sid}: {e}")
        traceback.print_exc()
