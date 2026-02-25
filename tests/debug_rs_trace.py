import sys, os
from rdkit import Chem
sys.path.append(os.path.abspath(os.path.dirname(__file__)))
from lib import sugar_sequence

def trace_signature(mol, ring_atoms):
    try:
        Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
        centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
        conf_map = {idx: conf for idx, conf in centers}
        print(f"Centers: {centers}")
        
        o_idx = None
        for idx in ring_atoms:
            if idx != -1 and mol.GetAtomWithIdx(idx).GetSymbol() == "O": o_idx=idx; break
        if o_idx is None: 
            print("No ring oxygen found in mapped ring!")
            return None
            
        o_atom = mol.GetAtomWithIdx(o_idx)
        neighbors = [n for n in o_atom.GetNeighbors() if n.GetIdx() in ring_atoms]
        if len(neighbors) != 2: 
            print(f"Oxygen has {len(neighbors)} neighbors in mapped ring!")
            return None
            
        n1, n2 = neighbors
        
        def score_c1(atom_idx):
            atom = mol.GetAtomWithIdx(atom_idx)
            s = 0
            for n in atom.GetNeighbors():
                if n.GetIdx() == o_idx: continue
                sym = n.GetSymbol()
                if sym in ('O', 'N', 'S', 'P'): s += 2
                if sym == 'C': s += 0.5
            return s
            
        s1 = score_c1(n1.GetIdx())
        s2 = score_c1(n2.GetIdx())
        print(f"n1={n1.GetIdx()} score={s1}, n2={n2.GetIdx()} score={s2}")
        
        c1_idx = None
        if s1 > s2: c1_idx = n1.GetIdx()
        elif s2 > s1: c1_idx = n2.GetIdx()
        else:
            def count_exo(atom_idx):
                atom = mol.GetAtomWithIdx(atom_idx)
                c = 0
                for n in atom.GetNeighbors():
                     if n.GetIdx() in ring_atoms: continue
                     if n.GetSymbol() != "H": c += 1
                return c
            e1 = count_exo(n1.GetIdx())
            e2 = count_exo(n2.GetIdx())
            if e1 >= e2: c1_idx = n1.GetIdx()
            else: c1_idx = n2.GetIdx()
            
        print(f"c1_idx selected: {c1_idx}")
        path = [c1_idx]
        curr = mol.GetAtomWithIdx(c1_idx)
        prev = o_idx
        
        for _ in range(len(ring_atoms)-2):
            found = False
            for n in curr.GetNeighbors():
                if n.GetIdx() in ring_atoms and n.GetIdx() != prev:
                    prev = curr.GetIdx()
                    curr = n
                    path.append(curr.GetIdx())
                    found = True
                    break
            if not found: 
                print(f"Broke early at {curr.GetIdx()}")
                break
                
        print(f"Path: {path}")
        sig = []
        for idx in path[1:]:
             sig.append(conf_map.get(idx, '?'))
        return tuple(sig)
    except Exception as e:
        print(f"Exception: {e}")
        return None

smi = "OC[C@H]1O[C@@H](O[C@@]2(CO)O[C@H](CO)[C@@H](O)[C@@H]2O)[C@H](O)[C@@H](O)[C@H]1O"
mol = Chem.MolFromSmiles(smi)
ri = mol.GetRingInfo()
ring = [list(r) for r in ri.AtomRings()][1] # fructose
print(f"Fructose target ring: {ring}")

clean_mol, mapping = sugar_sequence.isolate_sugar_ring(mol, ring)
mapped_ring = [mapping[idx] for idx in ring]
print(f"Mapped ring: {mapped_ring}")
for idx in mapped_ring:
    if idx != -1: print(f"  Atom {idx}: {clean_mol.GetAtomWithIdx(idx).GetSymbol()}")

sig = trace_signature(clean_mol, mapped_ring)
print(f"Trace sig: {sig}")

