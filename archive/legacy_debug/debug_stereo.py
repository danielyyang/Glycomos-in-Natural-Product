import sys, os
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np

sys.path.append(os.path.abspath(os.path.dirname(__file__)))
from lib import sugar_sequence

maltose = "O[C@H]1[C@H](O)[C@@H](O)[C@H](O[C@H]2[C@H](O)[C@@H](O)[C@H](O)[C@@H](CO)O2)[C@@H](CO)O1"
mol = Chem.MolFromSmiles(maltose)
ri = mol.GetRingInfo()
rings = [list(r) for r in ri.AtomRings()]
ring = rings[0] # Acceptor ring

clean_mol, mapping = sugar_sequence.isolate_sugar_ring(mol, ring)
mapped_ring = [mapping[idx] for idx in ring]

try:
    mol_calc = Chem.AddHs(clean_mol)
    AllChem.EmbedMolecule(mol_calc, randomSeed=42)
    conf = mol_calc.GetConformer()
    
    o_idx = None
    for idx in mapped_ring:
        if mol_calc.GetAtomWithIdx(idx).GetSymbol() == 'O':
            o_idx = idx
            break
            
    c1_idx = None
    for n_idx in [n.GetIdx() for n in mol_calc.GetAtomWithIdx(o_idx).GetNeighbors() if n.GetIdx() in mapped_ring]:
         n_atom = mol_calc.GetAtomWithIdx(n_idx)
         o_count = sum(1 for nn in n_atom.GetNeighbors() if nn.GetSymbol() == 'O')
         if o_count >= 2:
              c1_idx = n_idx
              break
              
    path = [c1_idx]
    curr = mol_calc.GetAtomWithIdx(c1_idx)
    prev = o_idx
    for _ in range(len(mapped_ring)-2):
        for n in curr.GetNeighbors():
            n_idx = n.GetIdx()
            if n_idx in mapped_ring and n_idx != prev:
                prev = curr.GetIdx()
                curr = n
                path.append(curr.GetIdx())
                break
                
    print(f"Path: {path}")

    coords = [conf.GetAtomPosition(idx) for idx in mapped_ring]
    center = np.mean(coords, axis=0)
    
    # Best fit plane normal using SVD
    points = np.array(coords) - center
    _, _, vh = np.linalg.svd(points)
    normal = vh[2, :] # normal vector
    
    print(f"Normal: {normal}")
    
    sig = []
    for c_idx in path[1:]:
        atom = mol_calc.GetAtomWithIdx(c_idx)
        subs = []
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() not in mapped_ring:
                subs.append(nbr)
                
        subs.sort(key=lambda x: x.GetAtomicNum(), reverse=True)
        sub_idx = subs[0].GetIdx()
        
        p_center = conf.GetAtomPosition(c_idx)
        p_sub = conf.GetAtomPosition(sub_idx)
        vec = np.array(p_sub) - np.array(p_center)
        
        projection = np.dot(vec, normal)
        direction = "U" if projection > 0 else "D"
        sig.append(direction)
        print(f"C{path.index(c_idx)+1} (Idx {c_idx}) -> Heavy Sub {mol_calc.GetAtomWithIdx(sub_idx).GetSymbol()} (Proj: {projection:.3f}) -> {direction}")
        
    print(f"Final Sig: {sig}")
    
except Exception as e:
    import traceback
    traceback.print_exc()

