# =============================================================================
# ⚠️ DEPRECATED — 此模块已被 sugar_utils.py 中的函数完全替代，保留仅供参考。
# ⚠️ DEPRECATED — This module has been fully superseded by functions in
# sugar_utils.py (find_nucleotides, extract_amino_acids_and_peptides).
# It is NOT used by main_pipeline.py. Retained for historical reference only.
# =============================================================================
from rdkit import Chem

# [EXPERT SMARTS DICTIONARY] 用于在脱糖骨架中检索核酸与肽链特征
SECONDARY_MOTIFS = {
    "Nucleotides": {
        # 嘌呤骨架 (Adenine, Guanine 等衍生物)
        "Purine_Base": "c1nc2c(n1)ncnc2",
        # 嘧啶骨架 (Cytosine, Uracil, Thymine 衍生物)
        "Pyrimidine_Base": "c1cncnc1",
        "Pyrimidine_Dione": "O=c1cc[nH]c(=O)[nH]1",
        # 磷酸基团
        "Phosphate": "P(=O)(-[OX2])"
    },
    "AminoAcids": {
        # 经典肽键 (Amide Bond) 连着手性碳
        "Peptide_Bond": "[NX3]C(=O)[CX4]",
        # 常见芳香族氨基酸残基特征 (如苯丙氨酸/酪氨酸残基)
        "Aromatic_Residue": "c1ccccc1CC(N)C(=O)"
    }
}

def extract_secondary_motifs(aglycan_smiles):
    """
    在脱糖基苷元 (Aglycan) 中扫描特定的次级代谢产物特征。
    """
    if not aglycan_smiles or aglycan_smiles == "NULL":
        return {"NUCLEOTIDES_SMILES": None, "AminoAcid_SMILES": None}
        
    mol = Chem.MolFromSmiles(aglycan_smiles)
    if not mol:
        return {"NUCLEOTIDES_SMILES": None, "AminoAcid_SMILES": None}

    found_nucleotides = []
    found_amino_acids = []

    # 检索核苷酸特征
    for name, smarts in SECONDARY_MOTIFS["Nucleotides"].items():
        pat = Chem.MolFromSmarts(smarts)
        if pat and mol.HasSubstructMatch(pat):
            found_nucleotides.append(name)

    # 检索氨基酸/肽链特征
    for name, smarts in SECONDARY_MOTIFS["AminoAcids"].items():
        pat = Chem.MolFromSmarts(smarts)
        if pat and mol.HasSubstructMatch(pat):
            found_amino_acids.append(name)

    return {
        "NUCLEOTIDES_SMILES": "|".join(found_nucleotides) if found_nucleotides else None,
        "AminoAcid_SMILES": "|".join(found_amino_acids) if found_amino_acids else None
    }

if __name__ == "__main__":
    # 测试 1: 尿苷二磷酸葡萄糖 (UDP-Glucose) 的 Aglycan (包含尿嘧啶和磷酸)
    # 假设糖链已被 Phase 2 切除，剩下的 Aglycan 带有 [15*]
    udp_aglycan = "[15*]OP(=O)(O)OP(=O)(O)OC[C@H]1O[C@H](O)[C@@H](O)[C@@H]1O"
    print("Testing UDP-Glucose Aglycan:", extract_secondary_motifs(udp_aglycan))
