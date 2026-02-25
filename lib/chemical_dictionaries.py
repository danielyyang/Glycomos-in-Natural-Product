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
    # D-Glucopyranose (Glc) - 基础赤道向构型
    "Glc": "[C:1]1([*:11])[C@H:2]([O:12])[C@@H:3]([O:13])[C@H:4]([O:14])[C@H:5]([CH2:6][O:16])[O:5]1",
    # D-Galactopyranose (Gal) - C4 差向异构（轴向）
    "Gal": "[C:1]1([*:11])[C@H:2]([O:12])[C@@H:3]([O:13])[C@@H:4]([O:14])[C@H:5]([CH2:6][O:16])[O:5]1",
    # D-Mannopyranose (Man) - C2 差向异构（轴向）
    "Man": "[C:1]1([*:11])[C@@H:2]([O:12])[C@@H:3]([O:13])[C@H:4]([O:14])[C@H:5]([CH2:6][O:16])[O:5]1",
    # L-Fucopyranose (Fuc) - 6-脱氧，L-半乳糖构型
    "Fuc": "[C:1]1([*:11])[C@@H:2]([O:12])[C@H:3]([O:13])[C@H:4]([O:14])[C@@H:5]([CH3:6])[O:5]1",
    # N-Acetyl-D-glucosamine (GlcNAc) - C2 被乙酰氨基取代
    "GlcNAc": "[C:1]1([*:11])[C@H:2]([N:12]C(=O)C)[C@@H:3]([O:13])[C@H:4]([O:14])[C@H:5]([CH2:6][O:16])[O:5]1",
    # D-Xylopyranose (Xyl) - 戊糖（无 C6）
    "Xyl": "[C:1]1([*:11])[C@H:2]([O:12])[C@@H:3]([O:13])[C@H:4]([O:14])[CH2:5][O:5]1",
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
# Format: [Target_Atom_In_Sugar:1]-[Modification_Atoms]>>[Target_Atom_In_Sugar:1][H]
_STRIPPING_SMIRKS_DEFS = {
    "Sulfated": "[O:1]S(=O)(=O)[O-]",  # O-Sulfate
    "Phosphated": "[O:1]P(=O)(O)O",    # O-Phosphate
    "Acetylated": "[O:1]C(=O)C",       # O-Acetyl
    "Methylated": "[O:1]C",            # O-Methyl (careful to run this after other C attachments if needed)
    "N-Acetylated": "[N:1]C(=O)C",     # N-Acetyl (e.g. GlcNAc -> GlcN)
    "N-Sulfated": "[N:1]S(=O)(=O)[O-]" # N-Sulfate
}

STRIPPING_REACTIONS = {}
for name, smirks in _STRIPPING_SMIRKS_DEFS.items():
    try:
        rxn = Chem.MolFromSmarts(smirks)
        if rxn:
            STRIPPING_REACTIONS[name] = rxn
    except Exception as e:
        print(f"Failed to load SMARTS for {name}: {e}")
