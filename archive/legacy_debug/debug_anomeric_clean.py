from rdkit import Chem
from rdkit.Chem import FindMolChiralCenters
import sys
sys.path.append('d:/Glycan_Database/lib')
import sugar_sequence

smiles = {
    "Pub_Mal": "C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@@H]2[C@H](O[C@H]([C@@H]([C@H]2O)O)CO)O)O)O)O)O",
    "User_Mal": "O[C@H]1[C@H](O)[C@@H](O)[C@H](O[C@H]2[C@H](O)[C@@H](O)[C@H](O)[C@@H](CO)O2)[C@@H](CO)O1",
    "Pub_Cel": "C([C@@H]1[C@H]([C@@H]([C@H]([C@@H](O1)O[C@@H]2[C@H](O[C@H]([C@@H]([C@H]2O)O)CO)O)O)O)O)O",
    "User_Cel": "O[C@H]1[C@H](O)[C@@H](O)[C@H](O[C@@H]2[C@H](O)[C@@H](O)[C@H](O)[C@@H](CO)O2)[C@@H](CO)O1"
}

for name, smi in smiles.items():
    print(f"--- {name} ---")
    mol = Chem.MolFromSmiles(smi)
    for ri, ring in enumerate(mol.GetRingInfo().AtomRings()):
        clean_mol, mapping = sugar_sequence.isolate_sugar_ring(mol, list(ring))
        if clean_mol is None: continue
        
        centers = dict(Chem.FindMolChiralCenters(clean_mol, includeUnassigned=True))
        mapped_ring = [mapping[idx] for idx in ring]
        
        # Get C1 and C5 the simple way:
        o_idx = next(idx for idx in mapped_ring if clean_mol.GetAtomWithIdx(idx).GetSymbol() == 'O')
        c1_idx = None
        for n_idx in [n.GetIdx() for n in clean_mol.GetAtomWithIdx(o_idx).GetNeighbors() if n.GetIdx() in mapped_ring]:
             n_atom = clean_mol.GetAtomWithIdx(n_idx)
             o_count = sum(1 for nn in n_atom.GetNeighbors() if nn.GetSymbol() == 'O')
             if o_count >= 2:
                  c1_idx = n_idx
                  break
        if c1_idx is None:
             c1_idx = [n.GetIdx() for n in clean_mol.GetAtomWithIdx(o_idx).GetNeighbors() if n.GetIdx() in mapped_ring][0]
        
        path = [c1_idx]
        curr = clean_mol.GetAtomWithIdx(c1_idx)
        prev = o_idx
        for _ in range(len(mapped_ring)-2):
            for n in curr.GetNeighbors():
                n_idx = n.GetIdx()
                if n_idx in mapped_ring and n_idx != prev:
                    prev = curr.GetIdx()
                    curr = n
                    path.append(curr.GetIdx())
                    break
                    
        c5_idx = path[-1]
        
        # Original Anomeric Logic
        c1_conf = centers.get(c1_idx, '?')
        c5_conf = centers.get(c5_idx, '?')
        
        is_beta = (c1_conf != c5_conf)
        anomer = 'b' if is_beta else 'a'
        print(f"Ring {ri}: C1 {c1_conf} | C5 {c5_conf} -> {anomer}")
