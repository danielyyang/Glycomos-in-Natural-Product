from rdkit import Chem
import sys

smis = {
    'Mal': 'O[C@H]1[C@H](O)[C@@H](O)[C@H](O[C@H]2[C@H](O)[C@@H](O)[C@H](O)[C@@H](CO)O2)[C@@H](CO)O1',
    'Cel': 'O[C@H]1[C@H](O)[C@@H](O)[C@H](O[C@@H]2[C@H](O)[C@@H](O)[C@H](O)[C@@H](CO)O2)[C@@H](CO)O1',
    'Suc': 'C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@]2([C@H]([C@@H]([C@H](O2)CO)O)O)CO)O)O)O)O',
    'Raf': 'C([C@@H]1[C@@H]([C@@H]([C@H]([C@H](O1)OC[C@@H]2[C@H]([C@@H]([C@H]([C@H](O2)O[C@]3([C@H]([C@@H]([C@H](O3)CO)O)O)CO)O)O)O)O)O)O)O',
    'Sta': 'C([C@@H]1[C@@H]([C@@H]([C@H]([C@H](O1)OC[C@@H]2[C@@H]([C@@H]([C@H]([C@H](O2)OC[C@@H]3[C@H]([C@@H]([C@H]([C@H](O3)O[C@]4([C@H]([C@@H]([C@H](O4)CO)O)O)CO)O)O)O)O)O)O)O)O)O)O'
}

for name, smi in smis.items():
    mol = Chem.MolFromSmiles(smi)
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    centers = dict(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
    print(f"--- {name} ---")
    
    for r in mol.GetRingInfo().AtomRings():
        # Get C1 and C5 inside the original ring
        o_idx = next(idx for idx in r if mol.GetAtomWithIdx(idx).GetSymbol() == 'O')
        c1_idx = None
        for n_idx in [n.GetIdx() for n in mol.GetAtomWithIdx(o_idx).GetNeighbors() if n.GetIdx() in r]:
             n_atom = mol.GetAtomWithIdx(n_idx)
             o_count = sum(1 for nn in n_atom.GetNeighbors() if nn.GetSymbol() == 'O')
             if o_count >= 2:
                  c1_idx = n_idx
                  break
        if c1_idx is None:
             c1_idx = [n.GetIdx() for n in mol.GetAtomWithIdx(o_idx).GetNeighbors() if n.GetIdx() in r][0]
        
        path = [c1_idx]
        curr = mol.GetAtomWithIdx(c1_idx)
        prev = o_idx
        for _ in range(len(r)-2):
            for n in curr.GetNeighbors():
                n_idx = n.GetIdx()
                if n_idx in r and n_idx != prev:
                    prev = curr.GetIdx()
                    curr = n
                    path.append(curr.GetIdx())
                    break
        c5_idx = path[-1]
        
        c1_conf = centers.get(c1_idx, '?')
        c5_conf = centers.get(c5_idx, '?')
        
        is_ketose = any(nbr.GetSymbol() == 'C' for nbr in mol.GetAtomWithIdx(c1_idx).GetNeighbors() if nbr.GetIdx() not in r)
        
        print(f"Ring -> C1: {c1_conf}, C5: {c5_conf}, Ketose: {is_ketose}")
