import sys
from rdkit import Chem

# Beta-D-Glucose
glc_smiles = "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O"
print(f"Glc (Ref): {glc_smiles}")

glc = Chem.MolFromSmiles(glc_smiles)

# Trace atoms to find C2 and C4
# SMILES: OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O
# O-C5-O-C1-C2-C3-C4-C5
# Indices in RDKit Mol (0-indexed) depends on SMILES parsing order.
# Let's map them.
# 0: O (of CH2OH?)
# 1: C (of CH2OH)
# 2: C (C5) - Ring
# 3: O (Ring O)
# 4: C (C1) - Ring
# 5: O (of C1-OH)
# 6: C (C2) - Ring
# 7: O (of C2-OH)
# 8: C (C3) - Ring
# 9: O (of C3-OH)
# 10: C (C4) - Ring
# 11: O (of C4-OH)

# Let's Verify connectivity
for atom in glc.GetAtoms():
    print(f"Atom {atom.GetIdx()} {atom.GetSymbol()} neighbors: {[x.GetIdx() for x in atom.GetNeighbors()]}")

# Based on "OC[C@H]1..."
# 0(O)-1(C)
# 1(C)-2(C) [C5]
# 2(C)-3(O) [Ring O]
# 3(O)-4(C) [C1]
# 4(C)-6(C) [C2]
# 6(C)-8(C) [C3]
# 8(C)-10(C) [C4]
# 10(C)-2(C) [Back to C5]

# So:
# C1 = 4
# C2 = 6
# C3 = 8
# C4 = 10
# C5 = 2

# Generate Mannose (C2 epimer)
man = Chem.Mol(glc)
man.GetAtomWithIdx(6).InvertChirality()
man_smiles = Chem.MolToSmiles(man, isomericSmiles=True)
print(f"Man (C2 inverted): {man_smiles}")

# Generate Galactose (C4 epimer)
gal = Chem.Mol(glc)
gal.GetAtomWithIdx(10).InvertChirality()
gal_smiles = Chem.MolToSmiles(gal, isomericSmiles=True)
print(f"Gal (C4 inverted): {gal_smiles}")
