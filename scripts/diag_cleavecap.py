"""
诊断 cleaveAndCap 碎片分类 bug (Diagnose fragment classification bug)
"""
import sys, os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "lib")))
from rdkit import Chem
from sugar_utils import find_mapped_sugar_units
from phase7_visualizer import identifySugarAtomZones

# Use a simple glycoside: rutin
smi = 'C=C1[C@@H]2CC[C@H]3[C@]4(C)C[C@H](O[C@H]5O[C@H](COC(=O)CC(C)C)[C@@H](OS(=O)(=O)O)[C@H](OS(=O)(=O)O)[C@H]5OC(=O)CC(C)C)CC(C(=O)O)(C(=O)O)[C@H]4CC[C@]3(C2)[C@H]1O'
mol = Chem.MolFromSmiles(smi)
units = find_mapped_sugar_units(mol)
print(f"Sugar units: {len(units)}")
for u in units:
    print(f"  {u['name']} ring_atoms={u['ring_atoms']}")

glycanAtoms, aglyconAtoms = identifySugarAtomZones(mol, units)
print(f"\nGlycanAtoms: {len(glycanAtoms)} -> {sorted(glycanAtoms)}")
print(f"AglyconAtoms: {len(aglyconAtoms)} -> {sorted(aglyconAtoms)}")

# Test FragmentOnBonds
boundaryBonds = []
for bond in mol.GetBonds():
    a1 = bond.GetBeginAtomIdx()
    a2 = bond.GetEndAtomIdx()
    if (a1 in glycanAtoms and a2 in aglyconAtoms) or (a1 in aglyconAtoms and a2 in glycanAtoms):
        boundaryBonds.append(bond.GetIdx())
print(f"\nBoundary bonds: {len(boundaryBonds)} -> {boundaryBonds}")

fragMol = Chem.FragmentOnBonds(mol, bondIndices=boundaryBonds, addDummies=True)
print(f"FragMol atoms: {fragMol.GetNumAtoms()} (original: {mol.GetNumAtoms()})")

# Check what GetMolFrags returns with fragsMolAtomMapping
fragsMolAtomMapping = []
fragMols = Chem.GetMolFrags(fragMol, asMols=True, sanitizeFrags=False, fragsMolAtomMapping=fragsMolAtomMapping)
print(f"\nFragments: {len(fragMols)}")
origNumAtoms = mol.GetNumAtoms()
for i, (frag, mapping) in enumerate(zip(fragMols, fragsMolAtomMapping)):
    numAtoms = frag.GetNumAtoms()
    # mapping[j] = atom index in fragMol for atom j in this fragment
    # Atoms with index < origNumAtoms are real atoms (original indices)
    # Atoms with index >= origNumAtoms are dummies inserted by FragmentOnBonds
    realOrigIndices = [m for m in mapping if m < origNumAtoms]
    dummyCount = sum(1 for m in mapping if m >= origNumAtoms)
    inGlycan = len(set(realOrigIndices) & glycanAtoms)
    inAglycon = len(set(realOrigIndices) & aglyconAtoms)
    hasSugarRing = False
    ri = frag.GetRingInfo()
    for ring in ri.AtomRings():
        if len(ring) in (5, 6):
            oCount = sum(1 for rIdx in ring if frag.GetAtomWithIdx(rIdx).GetAtomicNum() == 8)
            if oCount >= 1:
                hasSugarRing = True
    print(f"  Frag {i}: {numAtoms} atoms, {dummyCount} dummies, "
          f"inGlycan={inGlycan}, inAglycon={inAglycon}, sugarRing={hasSugarRing}")
    print(f"    mapping sample: {list(mapping[:5])}...")
    print(f"    SMILES: {Chem.MolToSmiles(frag)[:80]}")
