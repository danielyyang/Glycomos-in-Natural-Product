import traceback
from rdkit import Chem
import sys, os
sys.path.append(os.path.abspath(os.path.dirname(__file__)))
from lib import sugar_sequence

def log(msg):
    with open("d:/Glycan_Database/debug_log5.txt", "a", encoding="utf-8") as f:
        f.write(str(msg) + "\n")

if os.path.exists("d:/Glycan_Database/debug_log5.txt"):
    os.remove("d:/Glycan_Database/debug_log5.txt")

smi = "OC[C@H]1O[C@@H](O[C@@]2(CO)O[C@H](CO)[C@@H](O)[C@@H]2O)[C@H](O)[C@@H](O)[C@H]1O"
mol = Chem.MolFromSmiles(smi)
ri = mol.GetRingInfo()
rings = [list(r) for r in ri.AtomRings()]
ring = rings[1] # fructose
log(f"Fructose ring: {ring}")

try:
    log("Calling isolate_sugar_ring...")
    clean_mol, mapping = sugar_sequence.isolate_sugar_ring(mol, list(ring))
    log("Isolate returned.")
    if clean_mol and mapping:
        mapped_ring = [mapping[idx] for idx in ring]
        log(f"Mapped ring: {mapped_ring}")
        log("Calling get_rs_signature_core...")
        sig = sugar_sequence.get_rs_signature_core(clean_mol, mapped_ring)
        log(f"Signature: {sig}")
    else:
        log("No clean_mol or mapping")
        sig = sugar_sequence.get_rs_signature_core(mol, list(ring))
        
    log("Checking library...")
    base_name = "Hex" if len(ring)==6 else "Pen"
    log("Checking modifications...")
    mods = sugar_sequence.check_modifications(mol, ring)
    log(f"Mods: {mods}")
except Exception as e:
    log(f"CRASH: {e}")
    log(traceback.format_exc())
