import sys, os
from rdkit import Chem

def log(msg):
    with open("d:/Glycan_Database/debug_log3.txt", "a", encoding="utf-8") as f:
        f.write(msg + "\n")

if os.path.exists("d:/Glycan_Database/debug_log3.txt"):
    os.remove("d:/Glycan_Database/debug_log3.txt")

smi = "OC[C@H]1O[C@@H](O[C@@]2(CO)O[C@H](CO)[C@@H](O)[C@@H]2O)[C@H](O)[C@@H](O)[C@H]1O"
mol = Chem.MolFromSmiles(smi)

ri = mol.GetRingInfo()
rings = [list(r) for r in ri.AtomRings()]
ring = rings[1] # Fructose
ring_atoms = ring

log(f"Fructose ring: {ring}")

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
                    
bonds_to_cut = []
for exo_idx in list(exo_atoms):
    exo_atom = mol.GetAtomWithIdx(exo_idx)
    for nbr in exo_atom.GetNeighbors():
        nbr_idx = nbr.GetIdx()
        if nbr_idx not in ring_atoms and nbr_idx not in c6_idxs:
            bond = mol.GetBondBetweenAtoms(exo_idx, nbr_idx)
            if bond:
                bonds_to_cut.append(bond.GetIdx())
                
log(f"Bonds to cut: {bonds_to_cut}")

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
        log("Valid frag found.")
        log(Chem.MolToSmiles(frag))
        try:
            rwmol = Chem.RWMol(frag)
            for atom in rwmol.GetAtoms():
                if atom.GetAtomicNum() == 0:
                    nbrs = atom.GetNeighbors()
                    if nbrs:
                        nbr_sym = nbrs[0].GetSymbol()
                        log(f"Mutating dummy attached to {nbr_sym}")
                        if nbr_sym == 'C':
                            atom.SetAtomicNum(8)
                        elif nbr_sym in ('O', 'N', 'S'):
                            atom.SetAtomicNum(1)
                    atom.SetIsotope(0)
            log("Sanitizing...")
            Chem.SanitizeMol(rwmol)
            log("Success.")
        except Exception as e:
            log(f"Crash: {e}")

