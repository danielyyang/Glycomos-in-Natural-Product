import sys
from rdkit import Chem
from rdkit.Chem import FindMolChiralCenters

sugars = {
    "Glc": "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O", 
    "Gal": "OC[C@H]1O[C@@H](O)[C@@H](O)[C@@H](O)[C@@H]1O",
    "Man": "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O",
    "Fuc": "C[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O", 
    "Rha": "C[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O",
    "Xyl": "OC1CO[C@H](O)[C@H](O)[C@H]1O",
    "Ara": "C1[C@@H]([C@@H]([C@H]([C@H](O1)O)O)O)O", # L-Ara
    "GlcA": "O=C(O)[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O",
    "GalA": "O=C(O)[C@H]1O[C@@H](O)[C@@H](O)[C@@H](O)[C@@H]1O",
    "GlcNAc": "CC(=O)N[C@H]1[C@@H](O)O[C@H](CO)[C@@H](O)[C@@H]1O",
    "GalNAc": "CC(=O)N[C@H]1[C@@H](O)O[C@H](CO)[C@@H](O)[C@H]1O",
}

def get_rs_signature(mol):
    try:
        Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
        centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
        conf_map = {idx: conf for idx, conf in centers}
        
        ri = mol.GetRingInfo()
        if not ri.AtomRings(): return None
        ring = ri.AtomRings()[0]
        
        o_idx = None
        for idx in ring:
            if mol.GetAtomWithIdx(idx).GetSymbol() == "O": o_idx = idx; break
            
        o_atom = mol.GetAtomWithIdx(o_idx)
        neighbors = [n for n in o_atom.GetNeighbors() if n.GetIdx() in ring]
        n1, n2 = neighbors
        
        def is_anomeric(atom_idx):
            atom = mol.GetAtomWithIdx(atom_idx)
            for n in atom.GetNeighbors():
                if n.GetIdx() == o_idx: continue
                if n.GetSymbol() in ["O", "N"]: return True
            return False
            
        c1 = None; c5 = None
        if is_anomeric(n1.GetIdx()): c1 = n1; c5 = n2
        elif is_anomeric(n2.GetIdx()): c1 = n2; c5 = n1
        else:
             # Heuristic for Pentoses (C5 has H,H)
             # C1 has OH (O). C5 has H,H.
             def has_exo_heavy(atom_idx):
                atom = mol.GetAtomWithIdx(atom_idx)
                for n in atom.GetNeighbors():
                    if n.GetIdx() in ring: continue
                    if n.GetSymbol() != "H": return True
                return False
             if has_exo_heavy(n1.GetIdx()) and not has_exo_heavy(n2.GetIdx()): c1=n1; c5=n2
             elif has_exo_heavy(n2.GetIdx()) and not has_exo_heavy(n1.GetIdx()): c1=n2; c5=n1
             else: return None

        path = [c1.GetIdx()]
        curr = c1
        prev = o_idx
        # Walk ring
        for _ in range(len(ring)-2):
             found_next = False
             for n in curr.GetNeighbors():
                 if n.GetIdx() in ring and n.GetIdx() != prev:
                     prev = curr.GetIdx()
                     curr = n
                     path.append(curr.GetIdx())
                     found_next = True
                     break
             if not found_next: break
             
        # Extract sig
        sig = []
        # C2 is index 1
        for i in range(1, len(path)):
            sig.append(conf_map.get(path[i], '?'))
        return sig
    except:
        return None

print("RS_SIGS = {")
for name, s in sugars.items():
    mol = Chem.MolFromSmiles(s)
    if mol:
        sig = get_rs_signature(mol)
        print(f"    '{name}': {sig},")
print("}")
