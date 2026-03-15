"""
靶向验证测试: 验证 3 个 Bug 修复 (Targeted Verification Tests)
Test 1: CNP0426561.2 Name-Rescue fix
Test 2: CNP0207797.1 Sulfate + Acyl fix
"""
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "lib")))

from rdkit import Chem
from lib.sugar_utils import find_mapped_sugar_units
from lib.stereo_rescue import rescueViaName, replaceGenericLabelsInSequence

print("=" * 70)
print("  Test 1: CNP0426561.2 — Name-Rescue tokenizer fix")
print("=" * 70)

name = "6-C-Arabinopyranosyl-8-C-glucopyranosyltricin"
sequence = "Non ; Non"

# Step 1: 提取糖根词
inferred = rescueViaName(name, None)
print(f"  Name:     {name}")
print(f"  Seq:      {sequence}")
print(f"  Inferred: {inferred}")

# Step 2: 替换
rescued = replaceGenericLabelsInSequence(sequence, inferred)
print(f"  Rescued:  {rescued}")

# 验证
assert rescued != sequence, "FAIL: rescue should change the sequence"
assert "L-Ara" in rescued, f"FAIL: expected L-Ara in rescued seq, got: {rescued}"
assert "D-Glc" in rescued, f"FAIL: expected D-Glc in rescued seq, got: {rescued}"
print("  ✓ PASS: Name-Rescue correctly mapped both sugars\n")


print("=" * 70)
print("  Test 2: CNP0207797.1 — Sulfate + Acyl classification fix")
print("=" * 70)

smi = 'C=C1[C@@H]2CC[C@H]3[C@]4(C)C[C@H](O[C@H]5O[C@H](COC(=O)CC(C)C)[C@@H](OS(=O)(=O)O)[C@H](OS(=O)(=O)O)[C@H]5OC(=O)CC(C)C)CC(C(=O)O)(C(=O)O)[C@H]4CC[C@]3(C2)[C@H]1O'
mol = Chem.MolFromSmiles(smi)
units = find_mapped_sugar_units(mol)

print(f"  SMILES: {smi[:80]}...")
print(f"  Sugar units found: {len(units)}")

for u in units:
    print(f"  Name:  {u['name']}")
    print(f"  Mods:  {u['modifications']}")
    print(f"  Anomer: {u['anomeric_config']}")

    # 验证: 不应该包含 "Ac" (应该是 "Isovaleryl")
    assert "Ac" not in u['modifications'], \
        f"FAIL: 'Ac' should not be in mods, got: {u['modifications']}"
    # 验证: 应该包含 "Sulfate"
    assert "Sulfate" in u['modifications'], \
        f"FAIL: 'Sulfate' missing from mods, got: {u['modifications']}"
    # 验证: 应该包含 "Isovaleryl"
    assert "Isovaleryl" in u['modifications'], \
        f"FAIL: 'Isovaleryl' missing from mods, got: {u['modifications']}"

print("  ✓ PASS: Sulfate detected, Isovaleryl correctly classified (not Ac)\n")

print("=" * 70)
print("  Test 3: Edge cases — simple O-Ac should still work")
print("=" * 70)

# D-Glucose with O-Ac at C6
# 2,3,4-tri-O-acetyl-6-O-methyl-D-glucopyranose
acSmi = "CC(=O)O[C@@H]1[C@H](OC(C)=O)[C@@H](OC(C)=O)[C@H](O)[C@@H](COC)O1"
acMol = Chem.MolFromSmiles(acSmi)
if acMol:
    acUnits = find_mapped_sugar_units(acMol)
    for u in acUnits:
        print(f"  Name:  {u['name']}")
        print(f"  Mods:  {u['modifications']}")
        assert "O-Ac" in u['modifications'], \
            f"FAIL: simple O-Ac lost after fix! Got: {u['modifications']}"
    print("  ✓ PASS: simple O-Ac still detected correctly\n")
else:
    print("  SKIP: SMILES parse failed")

print("=" * 70)
print("  ALL TESTS PASSED")
print("=" * 70)
