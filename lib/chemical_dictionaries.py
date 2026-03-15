from rdkit import Chem

# ==============================================================================
# 1. Monosaccharide Library (Strict Stereochemical Whitelist)
# ==============================================================================
# We define these using generic SMILES. In SubstructMatch, RDKit allows substitution 
# on the OH groups (e.g., glycosidic bonds) by default when queried from SMILES,
# as long as we don't enforce exactly 1 degree or explicitly add Hs.
# We will enforce useChirality=True during matching.
# We must include both alpha and beta anomers, or use a non-stereo anomeric center 
# if we want to catch both anomers with a single query.
# Using 'C1(O)' without '@'/'@@' at C1 allows matching both alpha and beta anomers,
# while the rest of the ring stereochemistry precisely defines the sugar type (D-Glc, etc.).

_MONOSACCHARIDE_SMARTS_DEFS = {
    # === 基础己糖 (Standard Hexoses) ===
    # D-Glucopyranose (Glc) - 基础赤道向构型
    "Glc": "[C:1]1([*:11])[C@H:2]([O:12])[C@@H:3]([O:13])[C@H:4]([O:14])[C@H:5]([CH2:6][O:16])[O:5]1",
    # D-Galactopyranose (Gal) - C4 差向异构（轴向）
    "Gal": "[C:1]1([*:11])[C@H:2]([O:12])[C@@H:3]([O:13])[C@@H:4]([O:14])[C@H:5]([CH2:6][O:16])[O:5]1",
    # D-Mannopyranose (Man) - C2 差向异构（轴向）
    "Man": "[C:1]1([*:11])[C@@H:2]([O:12])[C@@H:3]([O:13])[C@H:4]([O:14])[C@H:5]([CH2:6][O:16])[O:5]1",

    # === 稀有己糖 (Rare Hexoses) ===
    # D-Allose (All) - C3 差向异构
    "All": "[C:1]1([*:11])[C@H:2]([O:12])[C@H:3]([O:13])[C@H:4]([O:14])[C@H:5]([CH2:6][O:16])[O:5]1",
    # D-Talose (Tal) - C2,C4 差向异构
    "Tal": "[C:1]1([*:11])[C@@H:2]([O:12])[C@H:3]([O:13])[C@@H:4]([O:14])[C@H:5]([CH2:6][O:16])[O:5]1",
    # D-Altrose (Alt) - C2,C3 差向异构
    "Alt": "[C:1]1([*:11])[C@@H:2]([O:12])[C@H:3]([O:13])[C@H:4]([O:14])[C@H:5]([CH2:6][O:16])[O:5]1",
    # D-Gulose (Gul) - C3,C4 差向异构
    "Gul": "[C:1]1([*:11])[C@H:2]([O:12])[C@H:3]([O:13])[C@@H:4]([O:14])[C@H:5]([CH2:6][O:16])[O:5]1",

    # === 脱氧糖 (Deoxy Sugars) ===
    # L-Fucopyranose (Fuc) - 6-脱氧，L-半乳糖构型
    "Fuc": "[C:1]1([*:11])[C@@H:2]([O:12])[C@H:3]([O:13])[C@H:4]([O:14])[C@@H:5]([CH3:6])[O:5]1",
    # L-Rhamnopyranose (Rha) - 6-脱氧-L-甘露糖构型
    "Rha": "[C:1]1([*:11])[C@H:2]([O:12])[C@H:3]([O:13])[C@@H:4]([O:14])[C@H:5]([CH3:6])[O:5]1",
    # D-Quinovose (Qui) - 6-deoxy-D-glucose
    "Qui": "[C:1]1([*:11])[C@H:2]([O:12])[C@@H:3]([O:13])[C@H:4]([O:14])[C@H:5]([CH3:6])[O:5]1",

    # === 氨基糖 (Amino Sugars) ===
    # N-Acetyl-D-glucosamine (GlcNAc) - C2 被乙酰氨基取代
    "GlcNAc": "[C:1]1([*:11])[C@H:2]([N:12]C(=O)C)[C@@H:3]([O:13])[C@H:4]([O:14])[C@H:5]([CH2:6][O:16])[O:5]1",
    # N-Acetyl-D-galactosamine (GalNAc) - C2 被乙酰氨基取代 (Gal 构型)
    "GalNAc": "[C:1]1([*:11])[C@H:2]([N:12]C(=O)C)[C@@H:3]([O:13])[C@@H:4]([O:14])[C@H:5]([CH2:6][O:16])[O:5]1",
    # N-Acetyl-D-mannosamine (ManNAc) - C2 被乙酰氨基取代 (Man 构型)
    "ManNAc": "[C:1]1([*:11])[C@@H:2]([N:12]C(=O)C)[C@@H:3]([O:13])[C@H:4]([O:14])[C@H:5]([CH2:6][O:16])[O:5]1",

    # === 糖醛酸 (Uronic Acids) — C6 = COOH ===
    # D-Glucuronic Acid (GlcA)
    "GlcA": "[C:1]1([*:11])[C@H:2]([O:12])[C@@H:3]([O:13])[C@H:4]([O:14])[C@H:5]([C:6](=O)[O:16])[O:5]1",
    # D-Galacturonic Acid (GalA)
    "GalA": "[C:1]1([*:11])[C@H:2]([O:12])[C@@H:3]([O:13])[C@@H:4]([O:14])[C@H:5]([C:6](=O)[O:16])[O:5]1",
    # L-Iduronic Acid (IdoA)
    "IdoA": "[C:1]1([*:11])[C@H:2]([O:12])[C@@H:3]([O:13])[C@H:4]([O:14])[C@@H:5]([C:6](=O)[O:16])[O:5]1",

    # === 戊糖 (Pentoses) ===
    # D-Xylopyranose (Xyl) - 戊糖（无 C6）
    "Xyl": "[C:1]1([*:11])[C@H:2]([O:12])[C@@H:3]([O:13])[C@H:4]([O:14])[CH2:5][O:5]1",
    # L-Arabinopyranose (Ara) - 皂苷/黄酮苷极常见
    "Ara": "[C:1]1([*:11])[C@@H:2]([O:12])[C@H:3]([O:13])[C@H:4]([O:14])[CH2:5][O:5]1",
    # D-Ribopyranose (Rib)
    "Rib": "[C:1]1([*:11])[C@H:2]([O:12])[C@H:3]([O:13])[C@H:4]([O:14])[CH2:5][O:5]1",
}

MONOSACCHARIDE_LIBRARY = {k: Chem.MolFromSmarts(v) for k, v in _MONOSACCHARIDE_SMARTS_DEFS.items()}

# Add wildcards to Anomeric centers for generic search fallback if needed
# But RDKit SMILES without explicit Hs normally match substituted versions inherently.

# ==============================================================================
# 2. Amino Acid Library (Strict Amide Bond Enforced)
# ==============================================================================
# We define SMARTS patterns for Amino Acids that STRICTLY require the amino/carboxylic 
# ends to be attached via an AMIDE bond pattern: [N;!R]-[C;!R](=O)
# If a free amine/carboxyl is detected at the termini of a standalone molecule, we allow it,
# but internally, it must peptide-bond.

# Legend:
# [NX3v3,NX4v4+] = Nitrogen (amine/amide)
# [CX4H] = Alpha Carbon (sp3, exactly 1 H for standard L-AA except Gly)
# [CX3](=O)[OX2,NX3] = Carboxyl or Amide Carbonyl

_AMINO_ACID_SMARTS_DEFS = {
    'Ala': '[NX3,NX4+][CX4H](C)[CX3](=O)[OX2,NX3]',
    'Arg': '[NX3,NX4+][CX4H](CCCNC(=[NH2,NH])[NH2,NH3+])[CX3](=O)[OX2,NX3]',
    'Asn': '[NX3,NX4+][CX4H](CC(=O)N)[CX3](=O)[OX2,NX3]',
    'Asp': '[NX3,NX4+][CX4H](CC(=O)[OH,O-])[CX3](=O)[OX2,NX3]',
    'Cys': '[NX3,NX4+][CX4H](CS)[CX3](=O)[OX2,NX3]',
    'Gln': '[NX3,NX4+][CX4H](CCC(=O)N)[CX3](=O)[OX2,NX3]',
    'Glu': '[NX3,NX4+][CX4H](CCC(=O)[OH,O-])[CX3](=O)[OX2,NX3]',
    'Gly': '[NX3,NX4+][CX4H2][CX3](=O)[OX2,NX3]',
    'His': '[NX3,NX4+][CX4H](Cc1cncn1)[CX3](=O)[OX2,NX3]',
    'Ile': '[NX3,NX4+][CX4H](C(C)CC)[CX3](=O)[OX2,NX3]',
    'Leu': '[NX3,NX4+][CX4H](CC(C)C)[CX3](=O)[OX2,NX3]',
    'Lys': '[NX3,NX4+][CX4H](CCCCN)[CX3](=O)[OX2,NX3]',
    'Met': '[NX3,NX4+][CX4H](CCSC)[CX3](=O)[OX2,NX3]',
    'Phe': '[NX3,NX4+][CX4H](Cc1ccccc1)[CX3](=O)[OX2,NX3]',
    'Pro': '[NX3,NX4+;r5]1[CX4H;r5](CCC1)[CX3](=O)[OX2,NX3]',
    'Ser': '[NX3,NX4+][CX4H](CO)[CX3](=O)[OX2,NX3]',
    'Thr': '[NX3,NX4+][CX4H](C(O)C)[CX3](=O)[OX2,NX3]',
    'Trp': '[NX3,NX4+][CX4H](Cc1c[nH]c2ccccc12)[CX3](=O)[OX2,NX3]',
    'Tyr': '[NX3,NX4+][CX4H](Cc1ccc(O)cc1)[CX3](=O)[OX2,NX3]',
    'Val': '[NX3,NX4+][CX4H](C(C)C)[CX3](=O)[OX2,NX3]',
}

AMINO_ACID_LIBRARY = {k: Chem.MolFromSmarts(v) for k, v in _AMINO_ACID_SMARTS_DEFS.items()}

# A strict sequence definition enforcing standard peptide backbone matching
STRICT_PEPTIDE_BOND_SMARTS = Chem.MolFromSmarts("[NX3][CX3](=O)")

# ==============================================================================
# 3. Nucleobase Library
# ==============================================================================
NUCLEOBASE_LIBRARY = {
    "Adenine": Chem.MolFromSmarts("n1cnc2c1ncnc2"),
    "Guanine": Chem.MolFromSmarts("n1cnc2c1nc(N)[nH]c2=O"),
    "Cytosine": Chem.MolFromSmarts("Nc1ccn([#6])c(=O)n1"), # Attached to sugar C
    "Thymine": Chem.MolFromSmarts("Cc1cn([#6])c(=O)[nH]c1=O"),
    "Uracil": Chem.MolFromSmarts("O=c1ccn([#6])c(=O)[nH]1")
}

# Phosphate rule for Nucleotides
PHOSPHATE_SMARTS = Chem.MolFromSmarts("P(=O)(O)(O)O")

# Reaction SMARTS (SMIRKS) for Modification Stripping
_STRIPPING_SMIRKS_DEFS = {
    # === 1. 无机酸修饰 (Inorganic Esters) ===
    "Sulfated": "[O:1]S(=O)(=O)[O-,OH]",        # 硫酸酯化 (O-Sulfate)
    "Phosphated": "[O:1]P(=O)([O-,OH])[O-,OH]", # 磷酸酯化 (O-Phosphate)
    # === 2. 基础与复杂脂肪酸酯 (Aliphatic Esters) ===
    "Acetylated": "[O:1]C(=O)[CH3]",            # 乙酰化 (O-Acetyl)
    "Formylated": "[O:1]C(=O)[H]",              # 甲酰化 (O-Formyl)
    "Malonylated": "[O:1]C(=O)CC(=O)[OH,O-]",   # 丙二酰化 (O-Malonyl)
    "Succinylated": "[O:1]C(=O)CCC(=O)[OH,O-]", # 琥珀酰化 (O-Succinyl)
    "Lactylated": "[O:1]C(C)C(=O)[OH,O-]",      # 乳酰化 (O-Lactyl)
    "Tigloylated": "[O:1]C(=O)C(=C)C",          # 巴豆酰/当归酰化 (O-Tigloyl/Angeloyl)
    # === 3. 芳香族高亮修饰 (Aromatic Esters - 天然产物重灾区) ===
    "Galloylated": "[O:1]C(=O)c1cc(O)c(O)c(O)c1",       # 没食子酰化 (O-Galloyl)
    "Benzoylated": "[O:1]C(=O)c1ccccc1",                # 苯甲酰化 (O-Benzoyl)
    "p-Coumaroylated": "[O:1]C(=O)/C=C/c1ccc(O)cc1",    # 对香豆酰化 (O-p-Coumaroyl)
    "Caffeoylated": "[O:1]C(=O)/C=C/c1ccc(O)c(O)c1",    # 咖啡酰化 (O-Caffeoyl)
    "Feruloylated": "[O:1]C(=O)/C=C/c1ccc(O)c(OC)c1",   # 阿魏酰化 (O-Feruloyl)
    "Sinapoylated": "[O:1]C(=O)/C=C/c1cc(OC)c(O)c(OC)c1", # 芥子酰化 (O-Sinapoyl)
    # === 4. 醚类修饰 (Ethers) ===
    "Methylated": "[O:1][CH3]",                 # 甲醚化 (O-Methyl, 已修复为安全末端甲基)
    # === 5. 氨基特种修饰 (N-Modifications) ===
    "N-Acetylated": "[N:1]C(=O)[CH3]",          # N-乙酰化 (N-Acetyl)
    "N-Glycolylated": "[N:1]C(=O)CO",           # N-羟乙酰化 (N-Glycolyl)
    "N-Sulfated": "[N:1]S(=O)(=O)[O-,OH]",      # N-硫酸酯化 (N-Sulfate)
    "N-Methylated": "[N:1][CH3]",               # N-甲基化 (N-Methyl)
    "N-Formylated": "[N:1]C(=O)[H]"             # N-甲酰化 (N-Formyl)
}

STRIPPING_REACTIONS = {}
for name, smirks in _STRIPPING_SMIRKS_DEFS.items():
    try:
        rxn = Chem.MolFromSmarts(smirks)
        if rxn:
            STRIPPING_REACTIONS[name] = rxn
    except Exception as e:
        print(f"Failed to load SMARTS for {name}: {e}")
