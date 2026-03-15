"""
Phase 2 碳丢失诊断: 用 Rutin 跑一遍当前切分逻辑
Diagnostic: run Rutin through current Phase 2 cleavage
"""
import sys, os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from rdkit import Chem
from lib.sugar_utils import find_mapped_sugar_units
from lib.cleavage_engine import cleave_glycan_aglycan

# Rutin: quercetin-3-O-[alpha-L-Rha-(1->6)-beta-D-Glc]
RUTIN_SMILES = "C[C@@H]1OC(OC[C@H]2OC(Oc3c(-c4ccc(O)c(O)c4)oc4cc(O)cc(O)c4c3=O)C(O)C(O)C2O)C(O)C(O)C1O"

mol = Chem.MolFromSmiles(RUTIN_SMILES)
if mol is None:
    # Try alternative Rutin SMILES
    RUTIN_SMILES = "O[C@@H]1[C@H](O)[C@@H](OC[C@H]2OC(Oc3c(-c4ccc(O)c(O)c4)oc5cc(O)cc(O)c5c3=O)[C@@H](O)[C@H](O)[C@@H]2O)O[C@@H](C)[C@@H]1O"
    mol = Chem.MolFromSmiles(RUTIN_SMILES)

print(f"Rutin SMILES: {RUTIN_SMILES}")
print(f"Mol valid: {mol is not None}")

if mol:
    totalC = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 6)
    totalAtoms = mol.GetNumAtoms()
    print(f"Total atoms: {totalAtoms}, Total C: {totalC}")

    # Phase 2: find sugar units
    units = find_mapped_sugar_units(mol)
    print(f"\nSugar units found: {len(units)}")
    for i, u in enumerate(units):
        print(f"  Unit {i}: ring_atoms={u['ring_atoms']}, name={u.get('name','?')}")

    # Phase 2: cleave
    glycan, aglycan = cleave_glycan_aglycan(mol, units)
    print(f"\nGlycan SMILES:  {glycan}")
    print(f"Aglycan SMILES: {aglycan}")

    # Carbon conservation check
    glycanC = 0
    aglycanC = 0
    if glycan and glycan != "NULL":
        gmol = Chem.MolFromSmiles(glycan)
        if gmol:
            glycanC = sum(1 for a in gmol.GetAtoms() if a.GetAtomicNum() == 6)
    if aglycan and aglycan != "NULL":
        amol = Chem.MolFromSmiles(aglycan)
        if amol:
            aglycanC = sum(1 for a in amol.GetAtoms() if a.GetAtomicNum() == 6)

    print(f"\nCarbon conservation:")
    print(f"  Original C: {totalC}")
    print(f"  Glycan C:   {glycanC}")
    print(f"  Aglycan C:  {aglycanC}")
    print(f"  Sum:        {glycanC + aglycanC}")
    print(f"  LOST:       {totalC - glycanC - aglycanC}")
    print(f"  Conservation OK: {totalC == glycanC + aglycanC}")
