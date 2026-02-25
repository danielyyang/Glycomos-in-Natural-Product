import traceback
from rdkit import Chem
import sys, os
sys.path.append(os.path.abspath(os.path.dirname(__file__)))
from lib import sugar_utils
from lib import sugar_sequence

def log(msg):
    with open("d:/Glycan_Database/debug_log.txt", "a", encoding="utf-8") as f:
        f.write(msg + "\n")

if os.path.exists("d:/Glycan_Database/debug_log.txt"):
    os.remove("d:/Glycan_Database/debug_log.txt")

log("Loading Sucrose...")
smi = "OC[C@H]1O[C@@H](O[C@@]2(CO)O[C@H](CO)[C@@H](O)[C@@H]2O)[C@H](O)[C@@H](O)[C@H]1O"
mol = Chem.MolFromSmiles(smi)

try:
    log("Running get_sugar_units...")
    units, atom_map = sugar_utils.get_sugar_units(mol)
    log("Success!")
except Exception as e:
    log(f"Crash in get_sugar_units: {e}")
    log(traceback.format_exc())
    
log("Done.")
