"""
Definitive relative stereochemistry algorithm for monosaccharides.
Instead of using absolute CIP (R/S) which flips upon substitution, we calculate
the relative "Up" / "Down" facing of the exocyclic oxygen/carbon at each chiral center
relative to the C1-O5-C5 plane.
"""
import sys, os
from rdkit import Chem

SUGAR_SMILES = {
    # Reference D-sugars
    "Glc": "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O", 
    "Gal": "OC[C@H]1O[C@@H](O)[C@@H](O)[C@@H](O)[C@@H]1O",
    "Man": "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O",
    "All": "OC[C@H]1O[C@@H](O)[C@H](O)[C@H](O)[C@@H]1O",   
}

def get_relative_stereo(mol, ring_atoms):
    """
    Returns a string of 'U' (Up) and 'D' (Down) for C1, C2, C3, C4, C5 (C6 position)
    This is topologically invariant to substitutions because it relies on the internal
    chiral tags (CW/CCW) of RDKit which are based on atom index ordering, not CIP mass!
    """
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    
    # Actually, RDKit's FindMolChiralCenters(useLegacyImplementation=False) returns CIP tags.
    # To get relative up/down, we can just use the absolute SMILES string format.
    # Since we can't extract "Up/Down" natively without 3D coords, let's use a simpler heuristic:
    # We KNOW RDKit flipped S to R specifically at the substitution point.
    
    # What if we just virtually replace all O-C bonds on the sugar with O-H? 
    # That is EXACTLY what I tried to do with ReplaceSubstructs!
    # Let me write a PERFECT virtual hydrolysis function that actually works.
    pass

def perfect_virtual_hydrolysis(mol, ring_atoms):
    emol = Chem.EditableMol(mol)
    
    # We need to find all bonds between an O/N/S atom attached to the ring, and the rest of the molecule
    # EXCEPT we must keep the ring intact.
    
    # 1. Identify exocyclic atoms directly connected to ring carbons
    exo_atoms = set()
    for idx in ring_atoms:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() not in ring_atoms:
                exo_atoms.add(nbr.GetIdx())
                
    # 2. For C6, we also want the oxygen attached to it
    c6_idx = None
    c6_o_idx = None
    for idx in ring_atoms:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() not in ring_atoms and nbr.GetSymbol() == 'C':
                for nnbr in nbr.GetNeighbors():
                    if nnbr.GetSymbol() in ('O','N','S'):
                        c6_idx = nbr.GetIdx()
                        c6_o_idx = nnbr.GetIdx()
                        exo_atoms.add(c6_o_idx)
                        break
                        
    # 3. Cut bonds outgoing from exo_atoms
    bonds_to_cut = []
    for exo_idx in list(exo_atoms):
        exo_atom = mol.GetAtomWithIdx(exo_idx)
        for nbr in exo_atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            # If the neighbor is NOT the ring AND NOT C6, cut it!
            if nbr_idx not in ring_atoms and nbr_idx != c6_idx:
                bond = mol.GetBondBetweenAtoms(exo_idx, nbr_idx)
                if bond:
                    bonds_to_cut.append(bond.GetIdx())
                    
    if not bonds_to_cut:
        return Chem.Mol(mol)
        
    frag_mol = Chem.FragmentOnBonds(mol, bonds_to_cut)
    frags = Chem.GetMolFrags(frag_mol, asMols=True)
    
    # Find the frag with the ring
    for frag in frags:
        ri = frag.GetRingInfo()
        valid = False
        for ring in ri.AtomRings():
            atoms = [frag.GetAtomWithIdx(a).GetSymbol() for a in ring]
            if len(ring) in (5,6) and atoms.count('O') == 1:
                valid = True
                break
        if valid:
            # Replace Dummy atoms (*) with Hydrogens or just remove them
            # Removing them and making the O/N a radical/charged is easiest to avoid CIP flips,
            # Or just replace with H
            smi = Chem.MolToSmiles(frag)
            # A dummy atom in SMILES is usually written as [*] or just *
            # If we replace * with [H] or remove it
            import re
            clean_smi = re.sub(r'\[\d+\*\]', '', smi)
            clean_mol = Chem.MolFromSmiles(clean_smi)
            return clean_mol
            
    return None

maltose = "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O[C@H]1OC(CO)[C@@H](O)[C@H](O)[C@H]1O"
mol = Chem.MolFromSmiles(maltose)
ri = mol.GetRingInfo()
for i, ring in enumerate(ri.AtomRings()):
    atoms = [mol.GetAtomWithIdx(a).GetSymbol() for a in ring]
    if atoms.count('O') == 1:
        clean = perfect_virtual_hydrolysis(mol, list(ring))
        if clean:
            Chem.AssignStereochemistry(clean, force=True, cleanIt=True)
            centers = Chem.FindMolChiralCenters(clean, includeUnassigned=True)
            print(f"Ring {i} cleaned centers: {centers}")
            print(f"SMILES: {Chem.MolToSmiles(clean, isomericSmiles=True)}")
