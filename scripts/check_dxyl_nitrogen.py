"""
在 debug sample 500 中搜索 D-Xyl + 含氮的化合物, 特别是糖环碳上有 N 的情况
Search debug sample 500 for D-Xyl with N on sugar ring carbons
"""
import sys, os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "lib")))

import pandas as pd
from rdkit import Chem
from lib import sugar_utils

df = pd.read_csv('reports/GlycoNP_Pipeline_Full_Cleaned.csv', low_memory=False, encoding='utf-8-sig')

# Search ALL D-Xyl rows for ones where N IS on the sugar ring carbon
xyl_rows = df[df['Sugar_Sequence'].astype(str).str.contains('D-Xyl', na=False)]
print(f"Total D-Xyl rows: {len(xyl_rows)}")

nitrogen_on_ring = []
checked = 0
for idx, row in xyl_rows.iterrows():
    smi = str(row.get('canonical_smiles', ''))
    if 'N' not in smi or smi in ('nan', ''):
        continue
    mol = Chem.MolFromSmiles(smi)
    if not mol:
        continue
    
    units, _ = sugar_utils.get_sugar_units(mol)
    for u in units:
        if 'Xyl' not in u.get('name', ''):
            continue
        ring = u['ring_atoms']
        ring_set = set(ring)
        for rIdx in ring:
            atom = mol.GetAtomWithIdx(rIdx)
            if atom.GetSymbol() != 'C':
                continue
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring_set:
                    continue
                if nbr.GetSymbol() == 'N':
                    cid = str(row.get('coconut_id', idx))
                    nitrogen_on_ring.append({
                        'id': cid,
                        'seq': str(row.get('Sugar_Sequence', '')),
                        'mods': str(row.get('Glycan_Modifications', '')),
                        'smi': smi[:100],
                        'n_idx': nbr.GetIdx(),
                        'n_hs': nbr.GetTotalNumHs(),
                        'n_degree': nbr.GetDegree()
                    })
    checked += 1
    if checked >= 5000:
        break

print(f"D-Xyl with N on sugar ring carbon (checked {checked}): {len(nitrogen_on_ring)}")
for item in nitrogen_on_ring[:20]:
    print(f"  ID={item['id']}  N_idx={item['n_idx']}  N_Hs={item['n_hs']}  N_deg={item['n_degree']}")
    print(f"    Seq={item['seq']}")
    print(f"    Mods={item['mods']}")
    print(f"    SMILES={item['smi']}")
    print()
