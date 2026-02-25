import sys
from rdkit import Chem
from rdkit.Chem import FindMolChiralCenters

# Copied from research
sugars = {
    # Aldohexoses (D-series mostly)
    "Glc": "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O", # beta-D-Glc (Ref)
    "Gal": "OC[C@H]1O[C@@H](O)[C@@H](O)[C@@H](O)[C@@H]1O", # beta-D-Gal
    "Man": "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O",   # beta-D-Man
    "All": "OC[C@H]1O[C@@H](O)[C@H](O)[C@H](O)[C@@H]1O",   # beta-D-Allose (C3 epimer of Glc? Verify)
           # Glc: C2(S), C3(R), C4(S), C5(R) [based on ring walk?]
           # All: C3 inverted vs Glc?
    "Alt": "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O", # Wait, Man is C2 epimer. 
           # Let's use the SMILES from the search result.
           # Search said: Alt: C([C@@H]1[C@H]([C@H]([C@@H]([C@@H](O1)O)O)O)O)O
           # I'll rely on the script to tell me the Sig.
    "Tal": "C([C@@H]1[C@@H]([C@@H]([C@@H]([C@@H](O1)O)O)O)O)O", # beta-D-Tal
    "Gul": "C([C@@H]1[C@H]([C@@H]([C@H]([C@@H](O1)O)O)O)O)O", # beta-D-Gul
    "Ido": "C([C@@H]1[C@H]([C@H]([C@@H]([C@@H](O1)O)O)O)O)O", # beta-D-Ido
    
    # Deoxy
    "Fuc": "C[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O", # alpha-L-Fuc
    "Rha": "C[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O",   # alpha-L-Rha (check C2 vs Fuc)
    "Qui": "C[C@@H]1[C@H]([C@@H]([C@H]([C@@H](O1)O)O)O)O", # beta-D-Qui (6-deoxy-Glc)
    
    # Pentoses
    "Xyl": "OC1CO[C@H](O)[C@H](O)[C@H]1O", # beta-D-Xyl
    "Ara": "C1[C@@H]([C@@H]([C@H]([C@H](O1)O)O)O)O", # beta-L-Ara
    "Rib": "C1[C@H]([C@H]([C@H]([C@@H](O1)O)O)O)O", # beta-D-Rib
    "Lyx": "C1[C@H]([C@@H]([C@@H]([C@@H](O1)O)O)O)O", # beta-D-Lyx
    
    # Acids
    "GlcA": "O=C(O)[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O", 
    "GalA": "O=C(O)[C@H]1O[C@@H](O)[C@@H](O)[C@@H](O)[C@@H]1O",
    "ManA": "O=C(O)[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O",
    "IdoA": "O=C(O)[C@H]1[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O", # L-Iduronic acid? (Isomer of GlcA/Ido?)
}

def get_rs_signature(mol):
    # Same logic as process_sugar.py but finding the single ring
    ri = mol.GetRingInfo()
    if not ri.AtomRings(): return None
    ring = ri.AtomRings()[0]
    
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    conf_map = {idx: conf for idx, conf in centers}
    
    o_idx = None
    for idx in ring:
        if mol.GetAtomWithIdx(idx).GetSymbol() == "O":
            o_idx = idx
            break
            
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
        # Fallback for Pentoses or if heuristic fails (Qui?)
        # For Pentoses, C5 is CH2 inside ring? No.
        # Pentose Pyranose: C1-C2-C3-C4-C5-O. No Exocyclic C on C5?
        # Xyl: C1..C5. Neighbors of C5 are O(ring), C4, H, H.
        # Neighbors of C1 are O(ring), C2, O(exo), H.
        # So C1 has Exo O. C5 does not.
        pass
        
    if c1 is None:
        # Use Exo atom check
        def has_exo_heavy(atom_idx):
            atom = mol.GetAtomWithIdx(atom_idx)
            for n in atom.GetNeighbors():
                if n.GetIdx() in ring: continue
                if n.GetSymbol() != "H": return True
            return False
        
        # C1 has Exo O. C5 has Exo C (Hex) or H (Xyl).
        # Both have exo heavy?
        # Xyl C1 has OH. C5 has H, H.
        pass

    # Basic walk relative to O
    # Let's just output the circular signature relative to O, then I can identify which is which manually?
    # Or try to enforce C1.
    
    # Just print the Chiral Centers in the ring
    ring_centers = []
    for idx in ring:
        if idx in conf_map:
            ring_centers.append((idx, conf_map[idx]))
            
    # Ordered walk
    path = [o_idx]
    curr = o_idx
    prev = None
    # Find n1 (start)
    curr = neighbors[0].GetIdx()
    path.append(curr)
    prev = o_idx
    
    while len(path) < len(ring) + 1:
        atom = mol.GetAtomWithIdx(curr)
        for n in atom.GetNeighbors():
            if n.GetIdx() in ring and n.GetIdx() != prev:
                prev = curr
                curr = n.GetIdx()
                path.append(curr)
                break
                
    # path is O -> C_a -> C_b ... -> C_z -> O
    # Collect configs
    sig = []
    for idx in path[1:-1]: # Skip O start and O end (overlap)
        sig.append(conf_map.get(idx, '?'))
        
    # We don't know direction (CW/CCW) or Start (C1 vs C5).
    # Return both directions.
    return sig, list(reversed(sig))

print("RS_LIBRARY = {")
for name, s in sugars.items():
    mol = Chem.MolFromSmiles(s)
    if not mol:
        print(f"  # Error parsing {name}")
        continue
    sig1, sig2 = get_rs_signature(mol)
    print(f"    '{name}': {sig1}, # Or {sig2}")
print("}")
