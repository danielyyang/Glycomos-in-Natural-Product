import traceback
from rdkit import Chem
import sys, os
sys.path.append(os.path.abspath(os.path.dirname(__file__)))
from lib import sugar_utils
from lib import sugar_sequence

def log(msg):
    with open("d:/Glycan_Database/debug_log4.txt", "a", encoding="utf-8") as f:
        f.write(str(msg) + "\n")

if os.path.exists("d:/Glycan_Database/debug_log4.txt"):
    os.remove("d:/Glycan_Database/debug_log4.txt")

smi = "OC[C@H]1O[C@@H](O[C@@]2(CO)O[C@H](CO)[C@@H](O)[C@@H]2O)[C@H](O)[C@@H](O)[C@H]1O"
mol = Chem.MolFromSmiles(smi)
ri = mol.GetRingInfo()
valid_rings = []
for ring in ri.AtomRings():
    atoms = [mol.GetAtomWithIdx(i) for i in ring]
    if [a.GetSymbol() for a in atoms].count("O") == 1:
        valid_rings.append(list(ring))

matched_units = []
Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
conf_map = {idx: conf for idx, conf in centers}

for sid, ring in enumerate(valid_rings):
    log(f"\nProcessing ring {sid}: {ring}")
    ring_o_idx = None
    for idx in ring:
        if mol.GetAtomWithIdx(idx).GetSymbol() == 'O':
            ring_o_idx = idx
            break
            
    o_atom = mol.GetAtomWithIdx(ring_o_idx)
    ring_neighbors = [n.GetIdx() for n in o_atom.GetNeighbors() if n.GetIdx() in ring]
    
    if len(ring_neighbors) != 2: 
        log("Not 2 neighbors")
        continue
    n1, n2 = ring_neighbors
    
    def score_c1(atom_idx):
        atom = mol.GetAtomWithIdx(atom_idx)
        score = 0
        for n in atom.GetNeighbors():
            idx = n.GetIdx()
            if idx in ring: continue
            sym = n.GetSymbol()
            if sym in ('O', 'N', 'S', 'P'): score += 2
            if sym == 'C': score += 0.5
        return score
        
    s1 = score_c1(n1)
    s2 = score_c1(n2)
    log(f"Scores: s1={s1}, s2={s2}")
    
    c1_idx = None
    if s1 > s2: c1_idx = n1
    elif s2 > s1: c1_idx = n2
    else:
        def count_exo(atom_idx):
            atom = mol.GetAtomWithIdx(atom_idx)
            return sum(1 for n in atom.GetNeighbors() if n.GetIdx() not in ring and n.GetSymbol() != 'H')
        if count_exo(n1) >= count_exo(n2): c1_idx = n1
        else: c1_idx = n2
        
    log(f"c1_idx selected: {c1_idx}")
    
    path = [c1_idx]
    curr = c1_idx
    prev = ring_o_idx
    log("Building path...")
    for _ in range(len(ring)-2):
        for n in mol.GetAtomWithIdx(curr).GetNeighbors():
            n_idx = n.GetIdx()
            if n_idx in ring and n_idx != prev:
                prev = curr
                curr = n_idx
                path.append(curr)
                break
                
    log(f"Path built: {path}")
    pos_map = {idx: i+1 for i, idx in enumerate(path)}
    
    oxygen_map = {}
    c6_idx = None
    log("Building oxygen field...")
    for idx, pos in list(pos_map.items()):
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            if nbr_idx not in ring:
                sym = nbr.GetSymbol()
                if sym in ('O', 'N', 'S'):
                    oxygen_map[nbr_idx] = pos
                elif sym == 'C' and pos == len(ring) - 1:
                    c6_idx = nbr_idx
                    pos_map[c6_idx] = pos + 1
                    for nnbr in nbr.GetNeighbors():
                        if nnbr.GetSymbol() in ('O', 'N', 'S') and nnbr.GetIdx() != idx:
                            oxygen_map[nnbr.GetIdx()] = pos + 1
                            break
                            
    log("Oxygen field built.")
    anomeric_idx = c1_idx
    
    c1_conf = conf_map.get(anomeric_idx, "?")
    ref_idx = path[-1]
    ref_conf = conf_map.get(ref_idx, "?")
    anomer_config = "?"
    if c1_conf in ('R', 'S') and ref_conf in ('R', 'S'):
        if c1_conf == ref_conf: anomer_config = "β"
        else: anomer_config = "α"
        
    log("Calling seq.identify_monosaccharide_v2...")
    try:
        base_name = sugar_sequence.identify_monosaccharide_v2(mol, ring)
        log(f"identify returned {base_name}")
    except Exception as e:
        log(f"CRASH IN IDENTIFY: {e}")
        log(traceback.format_exc())

log("Done!")
