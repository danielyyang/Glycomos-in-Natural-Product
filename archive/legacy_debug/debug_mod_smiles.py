from rdkit import Chem

test_cases = {
    "Alpha-D-GlcNAc": "CC(=O)N[C@H]1[C@@H](O)O[C@H](CO)[C@@H](O)[C@@H]1O",
    "Alpha-L-Rhamnose": "C[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O",
    "Beta-D-Glucuronic Acid": "O=C(O)[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O",
    "Kanamycin A": "C1[C@H]([C@@H]([C@H]([C@@H]([C@H]1N)O[C@@H]2[C@@H]([C@H]([C@@H]([C@H](O2)CN)O)O)O)O)O[C@@H]3[C@@H]([C@H]([C@@H]([C@H](O3)CO)O)N)O)N"
}

def traverse(mol, start_idx, avoid_atoms):
    visited = set([start_idx])
    q = [start_idx]
    while q:
        curr = q.pop(0)
        for nbr in mol.GetAtomWithIdx(curr).GetNeighbors():
            idx = nbr.GetIdx()
            if idx not in avoid_atoms and idx not in visited:
                visited.add(idx)
                q.append(idx)
    return visited

import networkx as nx

for name, smiles in test_cases.items():
    mol = Chem.MolFromSmiles(smiles)
    ri = mol.GetRingInfo()
    best_ring = None
    for r in ri.AtomRings():
        if len(r) == 6:
            best_ring = r
            break
    if not best_ring:
        continue
    
    print(f"\n--- {name} ---")
    for r_idx in best_ring:
        r_atom = mol.GetAtomWithIdx(r_idx)
        for nbr in r_atom.GetNeighbors():
            if nbr.GetIdx() not in best_ring:
                branch = traverse(mol, nbr.GetIdx(), set(best_ring))
                sub_smiles = Chem.MolFragmentToSmiles(mol, atomsToUse=list(branch))
                print(f"Branch at ring atom {r_idx} -> Exocyclic {nbr.GetSymbol()}: {sub_smiles}")
