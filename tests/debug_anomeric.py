from rdkit import Chem
from rdkit.Chem import FindMolChiralCenters
import sys
sys.path.append('d:/Glycan_Database')
from lib import sugar_utils

smiles = {
    "Pub_Mal": "C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@@H]2[C@H](O[C@H]([C@@H]([C@H]2O)O)CO)O)O)O)O)O",
    "User_Mal": "O[C@H]1[C@H](O)[C@@H](O)[C@H](O[C@H]2[C@H](O)[C@@H](O)[C@H](O)[C@@H](CO)O2)[C@@H](CO)O1",
    "Pub_Cel": "C([C@@H]1[C@H]([C@@H]([C@H]([C@@H](O1)O[C@@H]2[C@H](O[C@H]([C@@H]([C@H]2O)O)CO)O)O)O)O)O",
    "User_Cel": "O[C@H]1[C@H](O)[C@@H](O)[C@H](O[C@@H]2[C@H](O)[C@@H](O)[C@H](O)[C@@H](CO)O2)[C@@H](CO)O1"
}

for name, smi in smiles.items():
    print(f"--- {name} ---")
    mol = Chem.MolFromSmiles(smi)
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    centers = dict(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
    
    units, _ = sugar_utils.get_sugar_units(mol)
    for u in units:
        c1 = u['anomeric_idx']
        c5 = list(u['position_map'].keys())[list(u['position_map'].values()).index(5)]
        print(f"Ring {u['id']} ({u['name']}): C1 {centers.get(c1, '?')} | C5 {centers.get(c5, '?')} -> {u['anomeric_config']}")
