import sys, os
sys.path.append(os.path.abspath(os.path.dirname(__file__)))
from rdkit import Chem
from lib import sugar_sequence

def log(msg):
    with open("d:/Glycan_Database/debug_log2.txt", "a", encoding="utf-8") as f:
        f.write(msg + "\n")

if os.path.exists("d:/Glycan_Database/debug_log2.txt"):
    os.remove("d:/Glycan_Database/debug_log2.txt")

log("Loading Sucrose...")
smi = "OC[C@H]1O[C@@H](O[C@@]2(CO)O[C@H](CO)[C@@H](O)[C@@H]2O)[C@H](O)[C@@H](O)[C@H]1O"
mol = Chem.MolFromSmiles(smi)

ri = mol.GetRingInfo()
log("Finding rings...")
for sid, ring in enumerate(ri.AtomRings()):
    atoms = [mol.GetAtomWithIdx(i) for i in ring]
    if [a.GetSymbol() for a in atoms].count("O") == 1:
        log(f"Ring {sid}: {list(ring)}")
        try:
            log("Calling isolate...")
            ret = sugar_sequence.isolate_sugar_ring(mol, list(ring))
            log("Isolate success")
        except Exception as e:
            log(f"Crash in isolate: {e}")
log("Done")
