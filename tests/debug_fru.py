import sys, os
from rdkit import Chem
sys.path.append(os.path.abspath(os.path.dirname(__file__)))
from lib import sugar_sequence

print("Building Library...")
sugar_sequence.build_library()
print("Fru_Ref matched?", "Fru_Ref" in sugar_sequence.SUGAR_SMILES_LIB)

print("\nRS Library contents for Pen:")
for k, v in sugar_sequence.RS_LIBRARY.items():
    if len(k) == 3:
        print(f"  {k}: {v}")

smi = "OC[C@H]1O[C@@H](O[C@@]2(CO)O[C@H](CO)[C@@H](O)[C@@H]2O)[C@H](O)[C@@H](O)[C@H]1O"
mol = Chem.MolFromSmiles(smi)
ri = mol.GetRingInfo()
rings = [list(r) for r in ri.AtomRings()]
ring = rings[1] # fructose
print(f"Fructose ring: {ring}")

clean_mol, mapping = sugar_sequence.isolate_sugar_ring(mol, ring)
Chem.AssignStereochemistry(clean_mol, force=True, cleanIt=True)
mapped_ring = [mapping[idx] for idx in ring]
sig = sugar_sequence.get_rs_signature_core(clean_mol, mapped_ring)
print(f"Signature for Fructose in Sucrose: {sig}")

from rdkit.Chem import Draw
# We can also check C1 and linkage
from lib import sugar_utils
units, atom_map = sugar_utils.get_sugar_units(mol)
for u in units:
    print(f"Unit {u['id']}: name={u['name']}, anomeric_idx={u['anomeric_idx']}, ring_atoms={u['ring_atoms']}")
