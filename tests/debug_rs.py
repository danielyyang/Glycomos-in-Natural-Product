import sys
from rdkit import Chem
from rdkit.Chem import FindMolChiralCenters

ref_smiles = {
    "Glc_Beta": "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O",
    "Gal_Beta": "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O", 
    "Man_Beta": "OC[C@H]1O[C@@H](O)[C@@H](O)[C@@H](O)[C@@H]1O"
}

def get_ring_rs(mol):
    # Find ring atoms
    ri = mol.GetRingInfo()
    rings = ri.AtomRings()
    if not rings: return None
    ring = rings[0] # Assume 1 ring
    
    centers = FindMolChiralCenters(mol, includeUnassigned=True)
    # Map index -> config
    conf_map = {idx: conf for idx, conf in centers}
    
    # We need to order them C1-C2-C3-C4-C5.
    # Identify O first.
    o_idx = None
    for idx in ring:
        if mol.GetAtomWithIdx(idx).GetSymbol() == "O":
            o_idx = idx
            break
            
    if o_idx is None: return "No Ring O"
    
    # Identify C1 (neighbor of O with 2 O's attached or just Anomeric)
    # Simple: Neighbor of O that has another O neighbor (exocyclic)
    c1_idx = None
    c5_idx = None
    
    o_atom = mol.GetAtomWithIdx(o_idx)
    neighbors = o_atom.GetNeighbors()
    for n in neighbors:
        if n.GetIdx() not in ring: continue
        # Check if C1
        # C1 has another O neighbor usually (OH or OR)
        is_c1 = False
        for nn in n.GetNeighbors():
            if nn.GetIdx() == o_idx: continue
            if nn.GetSymbol() == "O":
                is_c1 = True
        
        if is_c1:
            c1_idx = n.GetIdx()
        else:
            c5_idx = n.GetIdx() # Typically C5 has CH2OH (C)
            
    # If both look like C1 (e.g. sucrose?), or ambiguity.
    # Glc/Gal/Man all have C1 with OH, C5 with CH2OH (C).
    # So C1 is connected to O (exocyclic), C5 connected to C (exocyclic).
    
    # Refine C1/C5 identification
    c1 = None
    c5 = None
    n1 = neighbors[0]
    n2 = neighbors[1]
    
    def is_anomeric(atom):
        # Has neighbor O that is not the ring O?
        for n in atom.GetNeighbors():
            if n.GetIdx() == o_idx: continue
            if n.GetSymbol() == "O": return True
        return False
        
    if is_anomeric(n1): c1=n1; c5=n2
    else: c1=n2; c5=n1
    
    # Walk C1 -> C2 -> C3 -> C4 -> C5
    # Traverse ring from C1 avoiding O
    
    current = c1
    prev = o_idx
    path = [c1.GetIdx()]
    
    for _ in range(4): # C2, C3, C4, C5
        for n in current.GetNeighbors():
            if n.GetIdx() in ring and n.GetIdx() != prev:
                prev = current.GetIdx()
                current = n
                path.append(current.GetIdx())
                break
                
    # Now get configs
    configs = []
    for idx in path:
        configs.append(conf_map.get(idx, '?'))
        
    return configs

for name, s in ref_smiles.items():
    mol = Chem.MolFromSmiles(s)
    print(f"{name}: {get_ring_rs(mol)}")
