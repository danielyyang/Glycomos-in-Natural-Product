from rdkit import Chem
import sys
sys.path.append('d:/Glycan_Database')
from lib import sugar_utils
from lib import sugar_sequence

def test_mol(name, smi):
    mol = Chem.MolFromSmiles(smi)
    units, _ = sugar_utils.get_sugar_units(mol)
    for u in units:
        print(f"[{name}] -> {u['name']} anomer: {u['anomeric_config']}")

smis = {
    'a-D-Glc': 'O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](CO)O1',
    'b-D-Glc': 'O[C@@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](CO)O1',
    'a-D-Gal': 'O[C@H]1[C@H](O)[C@@H](O)[C@@H](O)[C@@H](CO)O1',
    'b-D-Gal': 'O[C@@H]1[C@H](O)[C@@H](O)[C@@H](O)[C@@H](CO)O1',
    'a-D-Fru': 'O[C@]1(CO)O[C@H](CO)[C@@H](O)[C@@H]1O',
    'b-D-Fru': 'O[C@@]1(CO)O[C@H](CO)[C@@H](O)[C@@H]1O',
    'Pub_Sucrose': 'C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@]2([C@H]([C@@H]([C@H](O2)CO)O)O)CO)O)O)O)O',
    'Pub_Maltose': 'C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@@H]2[C@H](O[C@H]([C@@H]([C@H]2O)O)CO)O)O)O)O)O'
}

for name, smi in smis.items():
    test_mol(name, smi)
