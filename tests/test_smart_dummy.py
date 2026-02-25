"""
Test RWMol logic: convert dummy atoms to O or H based on neighbor.
"""
import sys, os
from rdkit import Chem

cellobiose = "O[C@H]1[C@H](O)[C@@H](O)[C@H](O[C@@H]2[C@H](O)[C@@H](O)[C@H](O)[C@@H](CO)O2)[C@@H](CO)O1"
mol = Chem.MolFromSmiles(cellobiose)

def fix_dummy_smart(frag_mol):
    rwmol = Chem.RWMol(frag_mol)
    for atom in rwmol.GetAtoms():
        if atom.GetAtomicNum() == 0:
            # find what it's attached to
            nbrs = atom.GetNeighbors()
            if not nbrs: continue
            nbr = nbrs[0]
            if nbr.GetSymbol() == 'C':
                atom.SetAtomicNum(8) # replace dummy with O
            elif nbr.GetSymbol() == 'O':
                atom.SetAtomicNum(1) # replace dummy with H
            atom.SetIsotope(0)
    # Important update implicit valences
    Chem.SanitizeMol(rwmol)
    return rwmol.GetMol()

# Fragment on all inter-ring bonds?
# Let's write the exact virtual hydrolysis bond finder:
ri = mol.GetRingInfo()
rings = [list(r) for r in ri.AtomRings()]
for i, ring in enumerate(rings):
    print(f"\n--- Ring {i} ---")
    
    exo_atoms = set()
    for idx in ring:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() not in ring:
                exo_atoms.add(nbr.GetIdx())
                
    c6_idx = None
    c6_o_idx = None
    for idx in ring:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() not in ring and nbr.GetSymbol() == 'C':
                for nnbr in nbr.GetNeighbors():
                    if nnbr.GetSymbol() in ('O','N','S'):
                        c6_idx = nbr.GetIdx()
                        c6_o_idx = nnbr.GetIdx()
                        exo_atoms.add(c6_o_idx)
                        break
                        
    bonds_to_cut = []
    for exo_idx in list(exo_atoms):
        exo_atom = mol.GetAtomWithIdx(exo_idx)
        for nbr in exo_atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            if nbr_idx not in ring and nbr_idx != c6_idx:
                bond = mol.GetBondBetweenAtoms(exo_idx, nbr_idx)
                if bond:
                    bonds_to_cut.append(bond.GetIdx())
                    
    frag_mol = Chem.FragmentOnBonds(mol, bonds_to_cut)
    frags = Chem.GetMolFrags(frag_mol, asMols=True)
    
    for frag in frags:
        f_ri = frag.GetRingInfo()
        valid = False
        for r in f_ri.AtomRings():
            atoms = [frag.GetAtomWithIdx(a).GetSymbol() for a in r]
            if len(r) in (5,6) and atoms.count('O') == 1:
                valid = True
                break
        if valid:
            print("Fragment before fix:", Chem.MolToSmiles(frag, isomericSmiles=True))
            clean = fix_dummy_smart(frag)
            Chem.AssignStereochemistry(clean, force=True, cleanIt=True)
            print("Fragment after fix: ", Chem.MolToSmiles(clean, isomericSmiles=True))
            centers = Chem.FindMolChiralCenters(clean)
            print(f"Centers: {centers}")
            # Isomeric match
            
