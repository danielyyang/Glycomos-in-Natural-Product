import sys
import os
sys.path.append(os.path.abspath(os.path.dirname(__file__)))
from rdkit import Chem
from process_sugar import generate_refined_sequence
from lib.sugar_utils import get_sugar_units, get_split_smiles

test_smiles = "CO[C@@H]1O[C@H](CO[C@@H]2O[C@H](COC(=O)OCC3c4ccccc4-c4ccccc43)[C@@H](OCc3ccccc3)[C@H](OCc3ccccc3)[C@H]2OC(=O)c2ccccc2)[C@@H](OCc2ccccc2)[C@H](OCc2ccccc2)[C@H]1OC(=O)c1ccccc1"
mol = Chem.MolFromSmiles(test_smiles)

print("1. Testing get_sugar_units")
try:
    units, a2s = get_sugar_units(mol)
    print(f"Found {len(units)} units")
except Exception as e:
    import traceback
    traceback.print_exc()

print("2. Testing generate_refined_sequence")
try:
    seq, mods = generate_refined_sequence(mol)
    print("SEQ:", seq)
    print("MODS:", mods)
except Exception as e:
    import traceback
    traceback.print_exc()

print("3. Testing get_split_smiles")
try:
    ag_sm, gl_sm = get_split_smiles(mol)
    print("AGLYCAN:", ag_sm)
    print("GLYCAN:", gl_sm)
except Exception as e:
    import traceback
    traceback.print_exc()

