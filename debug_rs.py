"""
Deep dive into the acceptor sugar of Maltose.
"""
import sys, os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__))))
from rdkit import Chem
from lib.sugar_sequence import isolate_sugar_ring, get_rs_signature_core, SUGAR_SMILES_LIB, RS_LIBRARY

maltose = "O[C@H]1[C@H](O)[C@@H](O)[C@H](O[C@H]2[C@H](O)[C@@H](O)[C@H](O)[C@@H](CO)O2)[C@@H](CO)O1"
mol = Chem.MolFromSmiles(maltose)
ri = mol.GetRingInfo()

print("Reference Glc:")
ref_glc = Chem.MolFromSmiles(SUGAR_SMILES_LIB["Glc"])
Chem.AssignStereochemistry(ref_glc, force=True, cleanIt=True)
glc_ring = ref_glc.GetRingInfo().AtomRings()[0]
glc_sig = get_rs_signature_core(ref_glc, list(glc_ring))
print(f"Centers: {Chem.FindMolChiralCenters(ref_glc)}")
print(f"Signature: {glc_sig}\n")

for i, ring in enumerate(ri.AtomRings()):
    atoms = [mol.GetAtomWithIdx(a).GetSymbol() for a in ring]
    if atoms.count('O') == 1:
        print(f"\n--- Ring {i} in Maltose ---")
        clean, mapping = isolate_sugar_ring(mol, list(ring))
        if clean and mapping:
            Chem.AssignStereochemistry(clean, force=True, cleanIt=True)
            mapped_ring = [mapping[idx] for idx in ring]
            sig = get_rs_signature_core(clean, mapped_ring)
            print(f"Centers: {Chem.FindMolChiralCenters(clean)}")
            print(f"SMILES: {Chem.MolToSmiles(clean, isomericSmiles=True)}")
            print(f"Signature: {sig}")
            if sig in RS_LIBRARY:
                print(f"  Matched: {RS_LIBRARY[sig]}")
            else:
                print("  NO MATCH")
                
# Add Fructose to SUGAR_SMILES_LIB if not present, and test Sucrose
SUGAR_SMILES_LIB["Fru"] = "OC[C@]1(O)O[C@H](CO)[C@@H](O)[C@@H]1O" # Beta-D-Fructofuranose
ref_fru = Chem.MolFromSmiles(SUGAR_SMILES_LIB["Fru"])
Chem.AssignStereochemistry(ref_fru, force=True, cleanIt=True)
fru_ring = ref_fru.GetRingInfo().AtomRings()[0]
fru_sig = get_rs_signature_core(ref_fru, list(fru_ring))
print(f"\nReference Fru:")
print(f"Centers: {Chem.FindMolChiralCenters(ref_fru)}")
print(f"Signature: {fru_sig}")
