from rdkit import Chem

amygdalin_smiles = "N#C[C@@H](c1ccccc1)O[C@@H]2O[C@H](CO[C@@H]3O[C@H](CO)[C@@H](O)[C@H](O)[C@H]3O)[C@@H](O)[C@H](O)[C@H]2O"

mol = Chem.MolFromSmiles(amygdalin_smiles)

glc_smarts = "[C:1]1([*:11])[C@H:2]([O:12])[C@@H:3]([O:13])[C@H:4]([O:14])[C@H:5]([CH2:6][O:16])[O:5]1"
glc_pat = Chem.MolFromSmarts(glc_smarts)

matches = mol.GetSubstructMatches(glc_pat, useChirality=True)
print(f"Number of matches with useChirality=True: {len(matches)}")
for m in matches:
    print(m)

matches_nochiral = mol.GetSubstructMatches(glc_pat, useChirality=False)
print(f"Number of matches with useChirality=False: {len(matches_nochiral)}")
for m in matches_nochiral:
    print(m)
