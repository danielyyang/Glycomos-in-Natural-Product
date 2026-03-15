import sys, os
from rdkit import Chem
sys.path.append(os.path.abspath(os.path.dirname(__file__)))
from lib import sugar_sequence

maltose = "O[C@H]1[C@H](O)[C@@H](O)[C@H](O[C@H]2[C@H](O)[C@@H](O)[C@H](O)[C@@H](CO)O2)[C@@H](CO)O1"
mol = Chem.MolFromSmiles(maltose)
ri = mol.GetRingInfo()
rings = [list(r) for r in ri.AtomRings()]
ring = rings[0] # Acceptor ring [2, 10, 8, 6, 4, 3]
print(f"Ring 0 (Acceptor): {ring}")

bonds_to_cut = []
ring_atoms = ring
exo_atoms = set()
for idx in ring_atoms:
    atom = mol.GetAtomWithIdx(idx)
    for nbr in atom.GetNeighbors():
        if nbr.GetIdx() not in ring_atoms:
            exo_atoms.add(nbr.GetIdx())
c6_idxs = []
for idx in ring_atoms:
    atom = mol.GetAtomWithIdx(idx)
    for nbr in atom.GetNeighbors():
        if nbr.GetIdx() not in ring_atoms and nbr.GetSymbol() == 'C':
            for nnbr in nbr.GetNeighbors():
                if nnbr.GetSymbol() in ('O','N','S'):
                    c6_idxs.append(nbr.GetIdx())
                    exo_atoms.add(nnbr.GetIdx())
                    break
for exo_idx in list(exo_atoms):
    exo_atom = mol.GetAtomWithIdx(exo_idx)
    for nbr in exo_atom.GetNeighbors():
        nbr_idx = nbr.GetIdx()
        if nbr_idx not in ring_atoms and nbr_idx not in c6_idxs:
            bond = mol.GetBondBetweenAtoms(exo_idx, nbr_idx)
            if bond: bonds_to_cut.append(bond.GetIdx())

frag_mol = Chem.FragmentOnBonds(mol, bonds_to_cut)
frags = Chem.GetMolFrags(frag_mol, asMols=True)

for i, frag in enumerate(frags):
    print(f"\nFrag {i}: {Chem.MolToSmiles(frag)}")
    matches = mol.GetSubstructMatches(frag)
    print(f"Matches count: {len(matches)}")
    for match in matches:
        print(f"  Match: {match}")
        overlap = sum(1 for t in match if t in ring_atoms)
        print(f"  Overlap count: {overlap}")
