"""Diagnose Glycan column rendering bug."""
import sys, os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "lib")))

from rdkit import Chem
from sugar_utils import find_mapped_sugar_units
from phase7_visualizer import identifySugarAtomZones, cleaveAndCap
from glycosidic_cleavage import cleaveWithConservation

TESTS = {
    "Rutin": "C[C@@H]1OC(OC[C@H]2OC(Oc3c(-c4ccc(O)c(O)c4)oc4cc(O)cc(O)c4c3=O)C(O)C(O)C2O)C(O)C(O)C1O",
    "GinsenosideRg1": "C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@]2(O)C[C@@H](C)O[C@@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@H]4O)C[C@H](O)[C@@H]5[C@@]3(C)CC[C@H](C5(C)C)O[C@@H]6O[C@H](CO)[C@@H](O)[C@H](O)[C@H]6O",
    # C-glycoside flavone (like image 1)
    "Vitexin_like": "OC1C(O)C(O)C(CO)OC1c1c(O)cc2oc(-c3ccc(O)cc3)cc(=O)c2c1O",
}

for name, smi in TESTS.items():
    mol = Chem.MolFromSmiles(smi)
    if not mol:
        print(f"{name}: Invalid SMILES"); continue
    units = find_mapped_sugar_units(mol)
    if not units:
        print(f"{name}: No sugars found"); continue
    sugarRingAtoms, substituentAtoms, aglyconAtomSet = identifySugarAtomZones(mol, units)
    glycanAtomSet = sugarRingAtoms | substituentAtoms
    fullH = mol.GetNumHeavyAtoms()

    print(f"\n{'='*60}")
    print(f"  {name}  (full={fullH} atoms)")
    print(f"{'='*60}")
    print(f"  zones: ring={len(sugarRingAtoms)} sub={len(substituentAtoms)} aglycon={len(aglyconAtomSet)} glycan={len(glycanAtomSet)}")

    # Direct from cleaveWithConservation
    try:
        gRaw, aRaw, meta = cleaveWithConservation(mol, units)
        gRawMol = Chem.MolFromSmiles(gRaw, sanitize=False) if gRaw and gRaw != "NULL" else None
        aRawMol = Chem.MolFromSmiles(aRaw, sanitize=False) if aRaw and aRaw != "NULL" else None
        gRawH = gRawMol.GetNumHeavyAtoms() if gRawMol else 0
        aRawH = aRawMol.GetNumHeavyAtoms() if aRawMol else 0
        print(f"  cleaveWithConservation:")
        print(f"    glycanSmi ({gRawH}H/{fullH}): {gRaw[:80] if gRaw else 'None'}")
        print(f"    aglyconSmi ({aRawH}H/{fullH}): {aRaw[:80] if aRaw else 'None'}")
        print(f"    bonds_cut: {meta.get('bonds_cut', 'N/A')}")
        if gRawH >= fullH * 0.85:
            print(f"    [!!] glycanSmi looks like FULL MOLECULE: {gRawH}/{fullH}")
    except Exception as e:
        print(f"  cleaveWithConservation ERROR: {e}")

    # From cleaveAndCap
    aFin, gFin = cleaveAndCap(mol, glycanAtomSet, aglyconAtomSet, minAglyconHeavyAtoms=3)
    gFinMol = Chem.MolFromSmiles(gFin, sanitize=False) if gFin else None
    aFinMol = Chem.MolFromSmiles(aFin, sanitize=False) if aFin else None
    gFinH = gFinMol.GetNumHeavyAtoms() if gFinMol else 0
    aFinH = aFinMol.GetNumHeavyAtoms() if aFinMol else 0
    print(f"  cleaveAndCap:")
    print(f"    aglyconFinal ({aFinH}H): {aFin[:80] if aFin else 'None'}")
    print(f"    glycanFinal ({gFinH}H): {gFin[:80] if gFin else 'None'}")

    if gFin and gFinH >= fullH * 0.85:
        print(f"    [BUG!] glycanFinal = FULL MOLECULE: {gFinH}/{fullH}")
    elif gFin:
        print(f"    [OK] glycan is {gFinH}/{fullH} ({gFinH/fullH*100:.0f}%% of full)")
    else:
        print(f"    [WARN] glycanFinal is None")
