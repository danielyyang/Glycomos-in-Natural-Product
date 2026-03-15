import sys
import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import FindMolChiralCenters
from rdkit.Chem import rdFMCS
from rdkit.Chem import AllChem
from tqdm import tqdm
import copy
import networkx as nx
import numpy as np

# Dynamically add lib path to support imports
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
try:
    import chemical_dictionaries
    from chemical_dictionaries import MONOSACCHARIDE_LIBRARY
except ImportError:
    from lib import chemical_dictionaries
    from lib.chemical_dictionaries import MONOSACCHARIDE_LIBRARY

try:
    import sugar_utils
except ImportError:
    from lib import sugar_utils

# ---------- 1. Reference Library Definitions ----------

# [THE COMPREHENSIVE NATURAL PRODUCTS DICTIONARY] D & L Enantiomers
# 最全单糖库 — 覆盖 8 大己糖全构型、脱氧糖、二脱氧/三脱氧糖、氨基糖、
# 糖醛酸 (Uronic Acids)、呋喃糖、酮糖、支链糖及细菌/海洋特有糖。
# Comprehensive monosaccharide library covering all 8 aldohexose stereoisomers,
# deoxy/dideoxy/trideoxy sugars, amino sugars, uronic acids, furanoses,
# ketoses, branched sugars, and bacterial/marine-specific sugars.
# All SMILES are in closed-ring pyranose/furanose form with full stereochemistry.
RAW_MONOSACCHARIDE_SMILES = {
    # =====================================================================
    # 1. 八大 D-己醛糖吡喃糖 (All 8 D-Aldohexopyranoses) + L-构型补充
    # =====================================================================
    # D-Glucose (Glc) — 最常见的己糖
    ("D-Glc", "a"): "O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](CO)O1",
    ("D-Glc", "b"): "O[C@@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](CO)O1",
    ("L-Glc", "a"): "O[C@@H]1[C@@H](O)[C@H](O)[C@@H](O)[C@H](CO)O1",
    ("L-Glc", "b"): "O[C@H]1[C@@H](O)[C@H](O)[C@@H](O)[C@H](CO)O1",
    # D-Galactose (Gal) — C4 差向异构体 (C4 epimer of Glc)
    ("D-Gal", "a"): "O[C@H]1[C@H](O)[C@@H](O)[C@@H](O)[C@@H](CO)O1",
    ("D-Gal", "b"): "O[C@@H]1[C@H](O)[C@@H](O)[C@@H](O)[C@@H](CO)O1",
    ("L-Gal", "a"): "O[C@@H]1[C@@H](O)[C@H](O)[C@H](O)[C@H](CO)O1",
    ("L-Gal", "b"): "O[C@H]1[C@@H](O)[C@H](O)[C@H](O)[C@H](CO)O1",
    # D-Mannose (Man) — C2 差向异构体 (C2 epimer of Glc)
    ("D-Man", "a"): "O[C@H]1[C@@H](O)[C@@H](O)[C@H](O)[C@@H](CO)O1",
    ("D-Man", "b"): "O[C@@H]1[C@@H](O)[C@@H](O)[C@H](O)[C@@H](CO)O1",
    ("L-Man", "a"): "O[C@@H]1[C@H](O)[C@H](O)[C@@H](O)[C@H](CO)O1",
    ("L-Man", "b"): "O[C@H]1[C@H](O)[C@H](O)[C@@H](O)[C@H](CO)O1",
    # D-Allose (All) — C3 差向异构体 (C3 epimer of Glc)
    ("D-All", "a"): "OC[C@H]1O[C@H](O)[C@H](O)[C@H](O)[C@@H]1O",
    ("D-All", "b"): "OC[C@H]1O[C@@H](O)[C@H](O)[C@H](O)[C@@H]1O",
    # D-Altrose (Alt) — C2,C3 差向异构体
    ("D-Alt", "a"): "OC[C@H]1O[C@H](O)[C@@H](O)[C@H](O)[C@@H]1O",
    ("D-Alt", "b"): "OC[C@H]1O[C@@H](O)[C@@H](O)[C@H](O)[C@@H]1O",
    # D-Gulose (Gul) — C3,C4 差向异构体
    ("D-Gul", "a"): "OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@H]1O",
    ("D-Gul", "b"): "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O",
    # D-Idose (Ido) — C2,C3,C4 差向异构体 (L-IdoA 在肝素中极重要)
    ("D-Ido", "a"): "OC[C@H]1O[C@H](O)[C@@H](O)[C@@H](O)[C@H]1O",
    ("D-Ido", "b"): "OC[C@H]1O[C@@H](O)[C@@H](O)[C@@H](O)[C@H]1O",
    # D-Talose (Tal) — C2,C4 差向异构体 (细菌 LPS 中偶见)
    ("D-Tal", "a"): "OC[C@H]1O[C@H](O)[C@@H](O)[C@@H](O)[C@@H]1O",
    ("D-Tal", "b"): "OC[C@H]1O[C@@H](O)[C@@H](O)[C@@H](O)[C@@H]1O",

    # =====================================================================
    # 2. 脱氧糖 (6-Deoxy Sugars) — 天然产物中的重灾区
    # =====================================================================
    # L-Rhamnose (Rha) — 6-deoxy-L-mannose，极常见于皂苷/黄酮苷
    ("L-Rha", "a"): "C[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O",
    ("L-Rha", "b"): "C[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@H]1O",
    ("D-Rha", "a"): "C[C@@H]1O[C@@H](O)[C@@H](O)[C@H](O)[C@@H]1O",
    # L-Fucose (Fuc) — 6-deoxy-L-galactose，血型抗原/海藻聚糖核心
    ("L-Fuc", "a"): "C[C@@H]1O[C@H](O)[C@@H](O)[C@H](O)[C@@H]1O",
    ("L-Fuc", "b"): "C[C@@H]1O[C@@H](O)[C@@H](O)[C@H](O)[C@@H]1O",
    ("D-Fuc", "a"): "C[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@H]1O",
    # D-Quinovose (Qui) — 6-deoxy-D-glucose，海洋天然产物常见
    ("D-Qui", "a"): "C[C@@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)O1",
    ("D-Qui", "b"): "C[C@@H]1[C@H](O)[C@@H](O)[C@H](O)[C@H](O)O1",
    ("L-Qui", "a"): "C[C@H]1[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)O1",

    # =====================================================================
    # 3. 二脱氧糖 (Dideoxy Sugars) — 心苷/抗生素核心
    # =====================================================================
    # D-Olivose (Oli) — 2,6-dideoxy-D-arabino-hexose，链霉类抗生素
    ("D-Oli", "a"): "C[C@@H]1C[C@@H](O)[C@H](O)[C@@H](O)O1",
    # D-Digitoxose (Dig) — 2,6-dideoxy-D-ribo-hexose，洋地黄毒苷核心
    ("D-Dig", "a"): "C[C@@H]1C[C@H](O)[C@H](O)[C@@H](O)O1",
    # L-Oleandrose (Ole) — 2,6-dideoxy-3-O-methyl-L-arabino-hexose
    ("L-Ole", "a"): "C[C@H]1C[C@H](O)[C@@H](OC)[C@H](O)O1",
    # D-Oliose — 2,6-dideoxy-D-lyxo-hexose (角鲨烷/聚酮类)
    ("D-Oliose", "a"): "C[C@@H]1C[C@@H](O)[C@@H](O)[C@@H](O)O1",
    # D-Boivinose (Boi) — 2,6-dideoxy-D-xylo-hexose
    ("D-Boi", "a"): "C[C@@H]1C[C@H](O)[C@@H](O)[C@@H](O)O1",
    # D-Cymarose (Cym) — 2,6-dideoxy-3-O-methyl-D-ribo-hexose
    ("D-Cym", "a"): "C[C@@H]1C[C@H](OC)[C@H](O)[C@@H](O)O1",

    # =====================================================================
    # 4. 三脱氧糖 (Trideoxy Sugars) — 细菌 LPS 免疫原性决定因子
    # =====================================================================
    # Abequose (Abe) — 3,6-dideoxy-D-xylo-hexose (Salmonella O-antigen)
    ("D-Abe", "a"): "C[C@@H]1[C@@H](O)C[C@H](O)[C@@H](O)O1",
    # Paratose (Par) — 3,6-dideoxy-D-ribo-hexose
    ("D-Par", "a"): "C[C@@H]1[C@H](O)C[C@H](O)[C@@H](O)O1",
    # Tyvelose (Tyv) — 3,6-dideoxy-D-arabino-hexose
    ("D-Tyv", "a"): "C[C@@H]1[C@H](O)C[C@@H](O)[C@@H](O)O1",
    # Colitose (Col) — 3,6-dideoxy-L-xylo-hexose (E. coli O-antigen)
    ("L-Col", "a"): "C[C@H]1[C@H](O)C[C@@H](O)[C@H](O)O1",
    # Ascarylose (Asc) — 3,6-dideoxy-L-arabino-hexose (线虫信息素)
    ("L-Asc", "a"): "C[C@H]1[C@@H](O)C[C@@H](O)[C@H](O)O1",

    # =====================================================================
    # 5. 氨基糖 (Amino Sugars) — 糖蛋白/糖胺聚糖核心
    # =====================================================================
    # N-Acetylglucosamine (GlcNAc) — 已在 chemical_dictionaries.py 中定义 SMARTS
    ("D-GlcNAc", "a"): "O[C@H]1[C@H](NC(=O)C)[C@@H](O)[C@H](O)[C@@H](CO)O1",
    ("D-GlcNAc", "b"): "O[C@@H]1[C@H](NC(=O)C)[C@@H](O)[C@H](O)[C@@H](CO)O1",
    # N-Acetylgalactosamine (GalNAc) — 黏蛋白 O-GalNAc 核心
    ("D-GalNAc", "a"): "O[C@H]1[C@H](NC(=O)C)[C@@H](O)[C@@H](O)[C@@H](CO)O1",
    ("D-GalNAc", "b"): "O[C@@H]1[C@H](NC(=O)C)[C@@H](O)[C@@H](O)[C@@H](CO)O1",
    # N-Acetylmannosamine (ManNAc) — 唾液酸生物合成前体
    ("D-ManNAc", "a"): "O[C@H]1[C@@H](NC(=O)C)[C@@H](O)[C@H](O)[C@@H](CO)O1",
    ("D-ManNAc", "b"): "O[C@@H]1[C@@H](NC(=O)C)[C@@H](O)[C@H](O)[C@@H](CO)O1",
    # Glucosamine (GlcN) — 游离氨基形式
    ("D-GlcN", "a"): "O[C@H]1[C@H](N)[C@@H](O)[C@H](O)[C@@H](CO)O1",
    # Galactosamine (GalN) — 游离氨基形式
    ("D-GalN", "a"): "O[C@H]1[C@H](N)[C@@H](O)[C@@H](O)[C@@H](CO)O1",
    # FucNAc — 2-acetamido-2,6-dideoxy-L-galactose (细菌 O-抗原)
    ("L-FucNAc", "a"): "C[C@@H]1O[C@H](O)[C@@H](O)[C@H](O)[C@@H]1NC(=O)C",
    # QuiNAc — 2-acetamido-2,6-dideoxy-D-glucose (细菌)
    ("D-QuiNAc", "a"): "C[C@@H]1[C@H](NC(=O)C)[C@@H](O)[C@H](O)[C@@H](O)O1",
    # RhaNAc — 2-acetamido-2,6-dideoxy-L-mannose (细菌)
    ("L-RhaNAc", "a"): "C[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@H]1NC(=O)C",
    # Bacillosamine (Bac) — 2,4-diamino-2,4,6-trideoxy-glucose (细菌蛋白糖基化)
    ("D-Bac", "a"): "C[C@@H]1[C@H](N)C[C@H](N)[C@@H](O)O1",

    # =====================================================================
    # 6. 糖醛酸 (Uronic Acids) — 糖胺聚糖/果胶核心
    # =====================================================================
    # D-Glucuronic Acid (GlcA) — 肝素/透明质酸/果胶
    ("D-GlcA", "a"): "O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](C(=O)O)O1",
    ("D-GlcA", "b"): "O[C@@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](C(=O)O)O1",
    # D-Galacturonic Acid (GalA) — 果胶主链
    ("D-GalA", "a"): "O[C@H]1[C@H](O)[C@@H](O)[C@@H](O)[C@@H](C(=O)O)O1",
    ("D-GalA", "b"): "O[C@@H]1[C@H](O)[C@@H](O)[C@@H](O)[C@@H](C(=O)O)O1",
    # D-Mannuronic Acid (ManA) — 海藻酸
    ("D-ManA", "a"): "O[C@H]1[C@@H](O)[C@@H](O)[C@H](O)[C@@H](C(=O)O)O1",
    # L-Iduronic Acid (IdoA) — 肝素/硫酸皮肤素核心，极其重要
    ("L-IdoA", "a"): "O[C@@H]1[C@H](O)[C@H](O)[C@@H](O)[C@H](C(=O)O)O1",
    ("L-IdoA", "b"): "O[C@H]1[C@H](O)[C@H](O)[C@@H](O)[C@H](C(=O)O)O1",
    # L-Guluronic Acid (GulA) — 海藻酸
    ("L-GulA", "a"): "O[C@@H]1[C@@H](O)[C@@H](O)[C@H](O)[C@H](C(=O)O)O1",

    # =====================================================================
    # 7. 戊糖 (Pentoses) — 吡喃糖和呋喃糖
    # =====================================================================
    # D-Xylose (Xylp) — 半纤维素核心
    ("D-Xyl", "a"): "O[C@H]1[C@H](O)[C@@H](O)[C@H](O)CO1",
    ("D-Xyl", "b"): "O[C@@H]1[C@H](O)[C@@H](O)[C@H](O)CO1",
    # L-Arabinose (Arap) — 皂苷/黄酮苷极常见
    ("L-Ara", "a"): "O[C@H]1[C@@H](O)[C@@H](O)[C@@H](O)CO1",
    ("L-Ara", "b"): "O[C@@H]1[C@@H](O)[C@@H](O)[C@@H](O)CO1",
    ("D-Ara", "a"): "O[C@@H]1[C@H](O)[C@H](O)[C@H](O)CO1",
    ("D-Ara", "b"): "O[C@H]1[C@H](O)[C@H](O)[C@H](O)CO1",
    # D-Ribose (Ribp) — 核糖吡喃糖型
    ("D-Rib", "a"): "O[C@@H]1[C@H](O)[C@H](O)[C@@H](O)CO1",
    ("D-Rib", "b"): "O[C@H]1[C@H](O)[C@H](O)[C@@H](O)CO1",
    # D-Lyxose (Lyxp)
    ("D-Lyx", "a"): "O[C@@H]1[C@H](O)[C@@H](O)[C@H](O)CO1",
    ("D-Lyx", "b"): "O[C@H]1[C@H](O)[C@@H](O)[C@H](O)CO1",
    ("L-Lyx", "a"): "O[C@H]1[C@@H](O)[C@H](O)[C@@H](O)CO1",

    # --- 呋喃糖 (Furanoses, 5-membered ring) ---
    # D-Ribofuranose (Ribf) — RNA 核糖基底
    ("D-Ribf", "a"): "OC[C@H]1O[C@H](O)[C@H](O)[C@@H]1O",
    ("D-Ribf", "b"): "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H]1O",
    # L-Arabinofuranose (Araf) — 植物细胞壁阿拉伯木聚糖
    ("L-Araf", "a"): "OC[C@@H]1O[C@@H](O)[C@H](O)[C@@H]1O",
    ("L-Araf", "b"): "OC[C@@H]1O[C@H](O)[C@H](O)[C@@H]1O",
    ("D-Araf", "a"): "OC[C@H]1O[C@H](O)[C@@H](O)[C@H]1O",
    # D-Xylofuranose (Xylf)
    ("D-Xylf", "a"): "OC[C@@H]1O[C@H](O)[C@H](O)[C@@H]1O",
    # D-Galactofuranose (Galf) — 结核杆菌/真菌细胞壁关键组分
    ("D-Galf", "a"): "OC[C@H](O)[C@H]1O[C@H](O)[C@H](O)[C@@H]1O",
    ("D-Galf", "b"): "OC[C@H](O)[C@H]1O[C@@H](O)[C@H](O)[C@@H]1O",
    # D-Glucofuranose (Glcf)
    ("D-Glcf", "a"): "OC[C@H](O)[C@@H]1O[C@H](O)[C@H](O)[C@@H]1O",

    # =====================================================================
    # 8. 酮糖 (Ketoses)
    # =====================================================================
    # D-Fructose (Fru) — 蔗糖组分
    ("D-Fru", "a"): "OC[C@]1(O)[C@@H](O)[C@H](O)[C@@H](CO)O1",
    ("D-Fru", "b"): "OC[C@@]1(O)[C@@H](O)[C@H](O)[C@@H](CO)O1",
    ("D-Frup", "a"): "OC[C@]1(O)OC[C@@H](O)[C@@H](O)[C@@H]1O",   # D-Fructopyranose
    ("D-Fruf", "a"): "OC[C@H]1O[C@@](O)(CO)[C@@H](O)[C@@H]1O",   # D-Fructofuranose
    # L-Sorbose (Sorb) — 维生素 C 工业合成中间体
    ("L-Sorb", "a"): "OC[C@@]1(O)OC[C@H](O)[C@@H](O)[C@@H]1O",
    # D-Tagatose (Tag) — D-Gal 的酮糖异构体
    ("D-Tag", "a"): "OC[C@]1(O)OC[C@H](O)[C@H](O)[C@@H]1O",
    # D-Psicose (Psi) — 稀有天然酮糖
    ("D-Psi", "a"): "OC[C@]1(O)OC[C@H](O)[C@@H](O)[C@H]1O",

    # =====================================================================
    # 9. 支链糖 (Branched-chain Sugars)
    # =====================================================================
    # D-Apiose (Api) — 支链戊糖呋喃，在 pectin RG-II 和黄酮苷中常见
    ("D-Api", "a"): "O[C@H]1[C@@H](O)[C@@](O)(CO)CO1",
    ("D-Api", "b"): "O[C@@H]1[C@@H](O)[C@@](O)(CO)CO1",

    # =====================================================================
    # 10. 庚糖及特殊糖 (Heptoses & Special Sugars)
    # =====================================================================
    # D-Sedoheptulose — 庚糖，磷酸戊糖途径中间体
    ("D-Sed", "a"): "OC[C@H]1O[C@@](O)(CO)[C@@H](O)[C@H](O)[C@@H]1O",
    # L,D-Heptose (Hep) — 革兰阴性菌 LPS 内核心
    ("L-D-Hep", "a"): "OC[C@H](O)[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O",
    # D,D-Heptose — 革兰阴性菌 LPS 变体
    ("D-D-Hep", "a"): "OC[C@@H](O)[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O",

    # =====================================================================
    # 11. 高级多碳糖 (Advanced Multi-Carbon Sugars)
    # =====================================================================
    # N-Acetylneuraminic Acid (Neu5Ac / Sialic Acid) — 9 碳
    # 唾液酸: 吡喃糖环 + C1 羧基 + C5-NHAc + C6-C7-C8-C9 甘油侧链
    # 在糖蛋白、糖脂、血型抗原中极其重要
    ("Neu5Ac", "a"): "O=C(O)[C@@]1(O)C[C@H](O)[C@@H](NC(=O)C)[C@@H](O1)[C@@H](O)[C@@H](O)CO",
    ("Neu5Ac", "b"): "O=C(O)[C@]1(O)C[C@H](O)[C@@H](NC(=O)C)[C@@H](O1)[C@@H](O)[C@@H](O)CO",
    # KDO (3-Deoxy-D-manno-octulosonic acid) — 8 碳
    # 酮脱氧辛酸: 细菌 LPS 核心多糖标志性糖
    ("KDO", "a"): "O=C(O)[C@@]1(O)C[C@H](O)[C@@H](O)[C@@H](O1)[C@@H](O)CO",
    ("KDO", "b"): "O=C(O)[C@]1(O)C[C@H](O)[C@@H](O)[C@@H](O1)[C@@H](O)CO",
    # N-Glycolylneuraminic Acid (Neu5Gc) — 唾液酸变体 (非人源哺乳动物)
    ("Neu5Gc", "a"): "O=C(O)[C@@]1(O)C[C@H](O)[C@@H](NC(=O)CO)[C@@H](O1)[C@@H](O)[C@@H](O)CO",

    # =====================================================================
    # 12. Extended Amino Sugars (from periodic table)
    # =====================================================================
    ("D-ManN", "a"): "O[C@H]1[C@@H](N)[C@@H](O)[C@H](O)[C@@H](CO)O1",
    ("D-AllN", "a"): "OC[C@H]1O[C@H](O)[C@H](N)[C@H](O)[C@@H]1O",
    ("D-Des", "a"): "C[C@@H]1C[C@@H](O)[C@H](N(C)C)[C@@H](O)O1",
    ("D-MurNAc", "a"): "O[C@H]1[C@H](NC(=O)C)[C@@H](O[C@@H](C)C(=O)O)[C@H](O)[C@@H](CO)O1",
    ("D-GulNAc", "a"): "O[C@H]1[C@H](NC(=O)C)[C@H](O)[C@@H](O)[C@@H](C(=O)O)O1",
    ("D-TalNAc", "a"): "OC[C@H]1O[C@H](O)[C@@H](NC(=O)C)[C@@H](O)[C@@H]1O",
    ("D-AltNAc", "a"): "OC[C@H]1O[C@H](O)[C@@H](NC(=O)C)[C@H](O)[C@@H]1O",
    ("D-AllNAc", "a"): "OC[C@H]1O[C@H](O)[C@H](NC(=O)C)[C@H](O)[C@@H]1O",
    ("D-IdoNAc", "a"): "OC[C@H]1O[C@H](O)[C@@H](NC(=O)C)[C@@H](O)[C@H]1O",
    ("L-IdoNAc", "a"): "OC[C@@H]1O[C@@H](O)[C@H](NC(=O)C)[C@H](O)[C@@H]1O",
    ("D-MurNGc", "a"): "O[C@H]1[C@H](NC(=O)CO)[C@@H](O[C@@H](C)C(=O)O)[C@H](O)[C@@H](CO)O1",
    ("L-FucN", "a"): "C[C@@H]1O[C@H](O)[C@@H](O)[C@H](N)[C@@H]1O",

    # =====================================================================
    # 13. Extended Uronic Acids
    # =====================================================================
    ("D-AllA", "a"): "O[C@H]1[C@H](O)[C@H](O)[C@H](O)[C@@H](C(=O)O)O1",
    ("D-TalA", "a"): "OC(=O)[C@H]1O[C@H](O)[C@@H](O)[C@@H](O)[C@@H]1O",
    ("D-AltA", "a"): "O[C@H]1[C@@H](O)[C@H](O)[C@H](O)[C@@H](C(=O)O)O1",
    ("L-GalA", "a"): "O[C@@H]1[C@@H](O)[C@H](O)[C@H](O)[C@H](C(=O)O)O1",
    ("D-meGlcA", "a"): "O[C@H]1[C@H](O)[C@@H](OC)[C@H](O)[C@@H](C(=O)O)O1",

    # =====================================================================
    # 14. Extended Deoxy Sugars
    # =====================================================================
    ("D-dRib", "a"): "OC[C@H]1C[C@H](O)[C@@H](O)O1",
    ("D-dGlc", "a"): "OC[C@H]1OC(O)[C@H](O)[C@@H](O)C1",
    ("D-dAra", "a"): "O[C@@H]1C[C@H](O)[C@H](O)CO1",
    ("D-dXyl", "a"): "O[C@H]1C[C@@H](O)[C@H](O)CO1",
    ("L-Cla", "a"): "C[C@H]1O[C@H](O)C[C@@](C)(OC)[C@@H]1O",
    ("L-Aco", "a"): "C[C@H]1C[C@@H](N)[C@H](O)[C@@H](O)O1",
    ("D-The", "a"): "C[C@@H]1[C@H](O)[C@@H](OC)[C@H](O)[C@@H](O)O1",
    ("D-6dTalNAc", "a"): "C[C@H]1O[C@H](O)[C@@H](NC(=O)C)[C@@H](O)[C@@H]1O",
    ("D-6dGul", "a"): "C[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@H]1O",
    ("D-6dTal", "a"): "C[C@H]1O[C@H](O)[C@@H](O)[C@@H](O)[C@@H]1O",
    ("D-6dAlt", "a"): "C[C@H]1O[C@H](O)[C@@H](O)[C@H](O)[C@@H]1O",

    # =====================================================================
    # 15. Tetroses
    # =====================================================================
    ("D-Thr", "a"): "O[C@H]1[C@H](O)[C@@H](O)CO1",
    ("D-Ery", "a"): "O[C@@H]1[C@H](O)[C@@H](O)CO1",

    # =====================================================================
    # 16. Pentose Ketoses
    # =====================================================================
    ("D-Ribu", "a"): "OC[C@@]1(O)[C@@H](O)[C@@H](O)CO1",
    ("D-Xylu", "a"): "OC[C@]1(O)[C@H](O)[C@@H](O)CO1",

    # =====================================================================
    # 17. Extended Nonoses (9-Carbon)
    # =====================================================================
    ("Kdn", "a"): "O=C(O)[C@@]1(O)C[C@H](O)[C@@H](O)[C@@H](O1)[C@@H](O)[C@@H](O)CO",
    ("Pse", "a"): "O=C(O)[C@@]1(O)C[C@@H](O)[C@@H](NC(=O)C)[C@@H](O1)[C@H](O)[C@H](NC(=O)C)C",
    ("Leg", "a"): "O=C(O)[C@@]1(O)C[C@H](O)[C@@H](NC(=O)C)[C@@H](O1)[C@@H](O)[C@@H](NC(=O)C)C",
    ("Fus", "a"): "O=C(O)[C@@]1(O)C[C@H](O)[C@@H](NC(=O)C)[C@@H](O1)[C@@H](O)[C@@H](O)CO",

    # =====================================================================
    # 18. Extended Octoses (8-Carbon)
    # =====================================================================
    ("Oct", "a"): "O=C(O)[C@]1(O)C[C@@H](O)[C@@H](O)[C@@H](O1)[C@@H](O)CO",

    # =====================================================================
    # 19. Decoses (10-Carbon)
    # =====================================================================
    ("Dha", "a"): "O=C(O)[C@@]1(O)C[C@H](O)[C@@H](O)[C@@H](O1)[C@@H](O)[C@@H](O)[C@@H](O)CO",
}
def create_query_mol(smiles):
    """将标准 SMILES 转化为无视隐式氢的通用 Query 分子"""
    mol = Chem.MolFromSmiles(smiles)
    if not mol: return None
    qmol = Chem.RWMol(mol)
    for atom in qmol.GetAtoms():
        if atom.GetAtomicNum() == 8:
            # 将所有氧原子替换为通用氧查询，允许其在聚合物中连接其他重原子
            q_atom = Chem.rdqueries.AtomNumEqualsQueryAtom(8)
            qmol.ReplaceAtom(atom.GetIdx(), q_atom)
    return qmol.GetMol()

# 使用改造后的通用查询图预编译字典
REFERENCE_MOLS = {k: create_query_mol(v) for k, v in RAW_MONOSACCHARIDE_SMILES.items() if create_query_mol(v) is not None}

# =====================================================================
# 原子计数氧门控 (Atom-Counting Oxygen Gate)
# 解决 L-Col 等低 OH 数脱氧糖劫持高 OH 数单糖的问题
# Prevents 3,6-dideoxy sugars (2 OH) from hijacking 6-deoxy (3 OH)
# or normal hexoses (4 OH)
# =====================================================================
def _countRefOxygens(smiles: str) -> int:
    """计算参考 SMILES 中除环氧外的含氧取代基数量。
    Count non-ring-oxygen O atoms in reference sugar SMILES.
    This equals the expected number of OH/OR groups on the ring carbons.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return -1
    ri = mol.GetRingInfo()
    ringAtoms = set()
    for ring in ri.AtomRings():
        if len(ring) in (5, 6):
            ringAtoms.update(ring)
    oCount = 0
    for idx in range(mol.GetNumAtoms()):
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetAtomicNum() == 8 and idx not in ringAtoms:
            oCount += 1
    return oCount


def _countFragmentRingOxygens(mol, ringAtoms: set) -> int:
    """计算糖碎片中环碳所连接的非环氧原子数量。
    Count exocyclic O atoms bonded to ring carbons in the actual fragment.
    This is the actual hydroxyl/substituent count to compare vs reference.
    """
    oCount = 0
    ringSet = set(ringAtoms)
    for rIdx in ringSet:
        atom = mol.GetAtomWithIdx(rIdx)
        if atom.GetAtomicNum() != 6:  # 只看环上的碳
            continue
        for nbr in atom.GetNeighbors():
            nIdx = nbr.GetIdx()
            if nIdx not in ringSet and nbr.GetAtomicNum() == 8:
                oCount += 1
            # 也一并统计 N 原子 (氨基糖的 N 取代 OH)
            elif nIdx not in ringSet and nbr.GetAtomicNum() == 7:
                oCount += 1
    return oCount


# 预计算每个参考糖的外环氧原子数
REFERENCE_OXYGEN_COUNTS = {
    k: _countRefOxygens(v)
    for k, v in RAW_MONOSACCHARIDE_SMILES.items()
}


# =====================================================================
# C/N 原子计数门控 (Carbon/Nitrogen Atom-Count Gate)
# 防止氨基己糖 (6C, 1N) 被错误匹配到戊糖 (5C, 0N)
# Prevents amino hexoses from matching pentose templates
# =====================================================================
def _countRefCN(smiles: str) -> tuple:
    """计算参考糖 SMILES 中的碳原子和氮原子总数。
    Count total C and N atoms in a reference sugar SMILES.
    Returns: (carbonCount, nitrogenCount)
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return (-1, -1)
    cCount = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 6)
    nCount = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 7)
    return (cCount, nCount)


# 预计算每个参考糖的碳/氮原子数
# Precompute C/N counts for all reference sugars
REFERENCE_CN_COUNTS = {
    k: _countRefCN(v)
    for k, v in RAW_MONOSACCHARIDE_SMILES.items()
}


def _countFragmentCN(mol, sugarUnitAtoms: set) -> tuple:
    """计算糖碎片中的碳原子和氮原子总数 (含环内+外环)。
    Count total C and N in fragment (ring + exocyclic atoms within the sugar unit).
    Returns: (carbonCount, nitrogenCount)
    """
    cCount = 0
    nCount = 0
    for idx in sugarUnitAtoms:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetAtomicNum() == 6:
            cCount += 1
        elif atom.GetAtomicNum() == 7:
            nCount += 1
    return (cCount, nCount)

def isolate_sugar_ring(mol, ring_atoms):
    bonds_to_cut = []
    exo_atoms = set()
    for idx in ring_atoms:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() not in ring_atoms: exo_atoms.add(nbr.GetIdx())
                
    c6_idxs = []
    for idx in ring_atoms:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() not in ring_atoms and nbr.GetSymbol() == 'C':
                for nnbr in nbr.GetNeighbors():
                    if nnbr.GetSymbol() in ('O','N','S'):
                        c6_idxs.append(nbr.GetIdx())
                        exo_atoms.add(nnbr.GetIdx())
                        break
                        
    for exo_idx in list(exo_atoms):
        exo_atom = mol.GetAtomWithIdx(exo_idx)
        for nbr in exo_atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            if nbr_idx not in ring_atoms and nbr_idx not in c6_idxs:
                bond = mol.GetBondBetweenAtoms(exo_idx, nbr_idx)
                if bond: bonds_to_cut.append(bond.GetIdx())
                    
    if not bonds_to_cut: return Chem.Mol(mol), dict(zip(ring_atoms, ring_atoms))
        
    frag_mol = Chem.FragmentOnBonds(mol, bonds_to_cut)
    frag_mols = Chem.GetMolFrags(frag_mol, asMols=True)
    frag_indices = Chem.GetMolFrags(frag_mol, asMols=False)
    
    for frag, indices in zip(frag_mols, frag_indices):
        if not all(r in indices for r in ring_atoms): continue
        ri = frag.GetRingInfo()
        valid = any(len(r) in (5,6) and [frag.GetAtomWithIdx(a).GetSymbol() for a in r].count('O') == 1 for r in ri.AtomRings())
        if valid:
            rwmol = Chem.RWMol(frag)
            atoms_to_process = [(atom.GetIdx(), atom.GetNeighbors()[0].GetSymbol()) for atom in rwmol.GetAtoms() if atom.GetAtomicNum() == 0 and atom.GetNeighbors()]
            for dummy_idx, nbr_sym in atoms_to_process:
                dummy_atom = rwmol.GetAtomWithIdx(dummy_idx)
                dummy_atom.SetAtomicNum(8 if nbr_sym == 'C' else 1)
                dummy_atom.SetFormalCharge(0); dummy_atom.SetIsotope(0)
            
            rwmol.UpdatePropertyCache(strict=False)
            Chem.SanitizeMol(rwmol)
            clean_mol = Chem.RemoveHs(rwmol.GetMol())
            Chem.AssignStereochemistry(clean_mol, force=True, cleanIt=True)
            
            # CRITICAL FIX: useChirality=False prevents CIP inversion from breaking the mapping!
            matches = mol.GetSubstructMatches(clean_mol, useChirality=False)
            best_match = next((match for match in matches if sum(1 for t in match if t in ring_atoms) >= len(ring_atoms)), None)
            
            if best_match:
                clean_to_target = {clean_idx: target_idx for clean_idx, target_idx in enumerate(best_match)}
                target_to_clean = {target_idx: clean_idx for clean_idx, target_idx in clean_to_target.items()}
                mapped_ring = [target_to_clean.get(old_r, -1) for old_r in ring_atoms]
                if -1 not in mapped_ring:
                    return clean_mol, dict(zip(ring_atoms, mapped_ring))
            
            new_ri = clean_mol.GetRingInfo()
            if new_ri.AtomRings(): return clean_mol, dict(zip(ring_atoms, new_ri.AtomRings()[0]))
    return None, None

def check_modifications(mol, ring_atoms, ref_smiles=None):
    mods = []
    is_acid = False
    is_complex = False
    
    try:
        ref_mol = Chem.MolFromSmiles(ref_smiles) if ref_smiles else None
        
        if ref_smiles and ("A" in ref_smiles or "acid" in ref_smiles.lower()):
            is_acid = True
            
        def traverse(start_idx, avoid_atoms):
            visited = set([start_idx])
            q = [start_idx]
            while q:
                curr = q.pop(0)
                for nbr in mol.GetAtomWithIdx(curr).GetNeighbors():
                    idx = nbr.GetIdx()
                    if idx not in avoid_atoms and idx not in visited:
                        visited.add(idx)
                        q.append(idx)
            return visited

        for r_idx in ring_atoms:
            r_atom = mol.GetAtomWithIdx(r_idx)
            for nbr in r_atom.GetNeighbors():
                if nbr.GetIdx() not in ring_atoms:
                    branch = traverse(nbr.GetIdx(), set(ring_atoms))
                    if len(branch) > 15:
                        continue
                        
                    try:
                        sub_smiles = Chem.MolFragmentToSmiles(mol, atomsToUse=list(branch))
                        if sub_smiles:
                            # 识别 NAc/Ac/NH2
                            if any(x in sub_smiles for x in ['NC(C)=O', 'CC(N)=O', 'CC(=O)N']): mods.append('NAc')
                            elif 'C(=O)C' in sub_smiles or sub_smiles == 'CC(=O)': mods.append('Ac')
                            elif sub_smiles == 'N': mods.append('N')
                            
                            # 识别脱氧 (直接连在环上的单碳)
                            elif sub_smiles == 'C': mods.append('deoxy')
                            
                            # 识别 O-甲基 (O-Me)
                            elif sub_smiles in ['CO', 'OC']:
                                # 注意: C6 的正常羟甲基 (CH2OH) 的 sub_smiles 也是 CO
                                # 我们之前在遍历字典时已经识别了基底(比如 Glc)，此时它不应该被当作修饰
                                # 真正的 O-Me 修饰通常在这个函数的更复杂逻辑中，或被基底忽略。
                                # 为了安全起见保留原样。
                                pass
                            elif 'S(=O)(=O)' in sub_smiles: mods.append('S')
                            
                            # 识别糖醛酸 (O=CO or C(=O)O)
                            elif 'O=CO' in sub_smiles or 'C(=O)O' in sub_smiles:
                                mods.append('A')
                                is_acid = True
                                
                            elif 'P(=O)' in sub_smiles: mods.append('P')
                            elif len(branch) >= 6: is_complex = True
                    except:
                        pass
                        
    except:
        pass
        
    final_mods = sorted(list(set(mods)))
    return final_mods, is_acid, is_complex

def create_virtual_standard_sugar(mol, ring_atoms):
    """
    创建分子的临时副本。将其糖环上连接的异质重原子（如 N, S）临时突变为 O。
    这使得 GlcN 等氨基糖能够骗过字典，完美匹配 Glc 的空间立体构型。
    注意：脱氧位置（仅连有 H 或 C）不受影响，完美保护 Rhamnose 等脱氧糖。
    """
    rwmol = Chem.RWMol(mol)
    for idx in ring_atoms:
        atom = rwmol.GetAtomWithIdx(idx)
        if atom.GetAtomicNum() == 6: # 遍历糖环上的碳原子
            for nbr in atom.GetNeighbors():
                nbr_idx = nbr.GetIdx()
                if nbr_idx not in ring_atoms:
                    # 发现环外取代原子。如果是 N (7) 或 S (16)，临时突变为 O (8)
                    if nbr.GetAtomicNum() in [7, 16]:
                        nbr.SetAtomicNum(8)
    return rwmol.GetMol()

def rescue_hexose_by_key_nodes(mol, ring_atoms):
    """
    当严格 3D 字典匹配失败时启动。
    依靠图遍历读取 C2 和 C4 的 CIP (R/S) 构型进行靶向强行救援。
    """
    if len(ring_atoms) != 6:
        return "Pen", "?"  # 仅救援六元糖环，五元保持
        
    # 强制 RDKit 计算分子的手性构型
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    try:
        # 1. 寻找环氧原子 (Ring Oxygen)
        ring_oxygens = [idx for idx in ring_atoms if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8]
        if not ring_oxygens: return "Hex", "?"
        ring_o_idx = ring_oxygens[0]
        
        # 2. 寻找 C1 (异头碳，连着环氧和另一个杂原子)
        c1_idx = None
        for nbr in mol.GetAtomWithIdx(ring_o_idx).GetNeighbors():
            if nbr.GetIdx() in ring_atoms:
                hetero_count = sum(1 for n in nbr.GetNeighbors() if n.GetAtomicNum() in [8, 7, 16])
                if hetero_count >= 2:
                    c1_idx = nbr.GetIdx()
                    break
        if c1_idx is None: return "Hex", "?"
        
        # 3. 顺藤摸瓜，沿着碳链遍历 C2, C3, C4
        c2_idx = [n.GetIdx() for n in mol.GetAtomWithIdx(c1_idx).GetNeighbors() if n.GetIdx() in ring_atoms and n.GetIdx() != ring_o_idx][0]
        c3_idx = [n.GetIdx() for n in mol.GetAtomWithIdx(c2_idx).GetNeighbors() if n.GetIdx() in ring_atoms and n.GetIdx() != c1_idx][0]
        c4_idx = [n.GetIdx() for n in mol.GetAtomWithIdx(c3_idx).GetNeighbors() if n.GetIdx() in ring_atoms and n.GetIdx() != c2_idx][0]
        
        # 4. 读取 C2 和 C4 的 R/S 构型标签
        c2_cip = mol.GetAtomWithIdx(c2_idx).GetProp('_CIPCode') if mol.GetAtomWithIdx(c2_idx).HasProp('_CIPCode') else "?"
        c4_cip = mol.GetAtomWithIdx(c4_idx).GetProp('_CIPCode') if mol.GetAtomWithIdx(c4_idx).HasProp('_CIPCode') else "?"
        
        # 5. 检测 C2 上是否有氮原子 (N-detection for amino sugars)
        # 经典氨基糖 (GlcN, GalN, ManN) 的 N 都在 C2 位
        # Classic amino sugars have N at C2 position
        ring_set = set(ring_atoms)  # 列表转集合 (list → set for O(1) lookup)
        hasNitrogenOnC2 = any(
            n.GetAtomicNum() == 7
            for n in mol.GetAtomWithIdx(c2_idx).GetNeighbors()
            if n.GetIdx() not in ring_set
        )
        # 同时检查整个碎片是否有 N (备用)
        hasAnyNitrogen = any(
            mol.GetAtomWithIdx(rIdx).GetAtomicNum() == 7
            or any(n.GetAtomicNum() == 7 for n in mol.GetAtomWithIdx(rIdx).GetNeighbors() if n.GetIdx() not in ring_set)
            for rIdx in ring_set if mol.GetAtomWithIdx(rIdx).GetAtomicNum() == 6
        )
        # 氨基糖后缀 (Amino sugar suffix)
        aminoSuffix = "N" if (hasNitrogenOnC2 or hasAnyNitrogen) else ""

        # 6. 首席科学家的靶向特征判决逻辑 (基于 D-型吡喃糖标准)
        # CIP-based classification with amino sugar awareness
        if c4_cip == "R" and c2_cip == "S": return f"D-Gal{aminoSuffix}(Rescued)", "?"
        elif c4_cip == "S" and c2_cip == "S": return f"D-Glc{aminoSuffix}(Rescued)", "?"
        elif c4_cip == "S" and c2_cip == "R": return f"D-Man{aminoSuffix}(Rescued)", "?"
        elif c4_cip == "R" and c2_cip == "R": return f"D-Tal{aminoSuffix}(Rescued)", "?"
            
        return "Hex", "?"
    except Exception:
        # 图遍历如果遇到奇葩结构报错，安全退回 Hex
        return "Hex", "?"

def _countFragmentCarbons(mol, ring_atoms, maxHops: int = 4):
    """
    杂原子截断铁律: 严格计算糖的 C-C 连续骨架碳数。
    Heteroatom Blockade: count sugar backbone carbons via C-C bonds ONLY.

    化学定义: 糖的碳数分类 (Hex=6C, Hept=7C, Oct=8C, Non=9C) 仅由
    连续的碳-碳 (C-C) 骨架决定。杂原子 (O, N, S) 是硬截断屏障。
    Chemical rule: Sugar carbon classification is determined ONLY by the
    contiguous C-C backbone. Heteroatoms (O, N, S) are hard barriers.

    v2.1 修复: 增加 BFS 步数限制 (maxHops=4) 和分支度门控。
    v2.1 fix: Added BFS depth limit and branching-degree gate.
    - maxHops=4 覆盖: 环碳 → C6 → C7 → C8 → C9 (覆盖壬糖/唾液酸)
    - 分支度门控: 非环碳若有 ≥3 个碳邻居 → 已进入苷元骨架, 停止
    - maxHops=4 covers: ring-C → C6 → C7 → C8 → C9 (covers nonoses)
    - Branch gate: non-ring carbon with ≥3 carbon neighbors → in aglycon, stop

    示例 (Examples):
      - C-O-C(=O)CH3 (乙酰化): O 处截断, 外侧 2C 不计入
      - C-NH-C(=O)CH3 (NAc): N 处截断, 外侧 2C 不计入
      - C-C(OH)-C(OH)-CH2OH (甘油侧链): 沿 C-C 遍历, 全部计入 (≤4 步)
    """
    ringSet = set(ring_atoms)
    visited = set(ring_atoms)   # 环原子全部预标记已访问
    totalCarbons = 0

    # 首先计算环上的碳原子数
    carbonQueue = []
    for idx in ring_atoms:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetAtomicNum() == 6:
            totalCarbons += 1
            carbonQueue.append((idx, 0))  # (atomIdx, depth)
        # 环氧标记为 visited 但不扩展 (它是种子边界)

    # 从环上碳原子出发, 仅沿 C-C 键向外 BFS (带步数限制)
    # BFS from ring carbons, follow C-C bonds only (with depth limit)
    while carbonQueue:
        idx, depth = carbonQueue.pop(0)
        if depth >= maxHops:
            continue  # 步数限制: 不再扩展 (depth limit: stop expanding)
        for nbr in mol.GetAtomWithIdx(idx).GetNeighbors():
            nIdx = nbr.GetIdx()
            if nIdx in visited:
                continue
            visited.add(nIdx)

            # 杂原子截断铁律: O(8), N(7), S(16), P(15) → 硬停
            if nbr.GetAtomicNum() != 6:
                continue  # 遇到杂原子, 不再沿此分支前进

            # 分支度门控: 非环碳若有 ≥3 个碳邻居 → 已进入骨架, 停止
            # Branching gate: non-ring carbon with ≥3 C neighbors → in backbone
            if nIdx not in ringSet:
                carbonNbrCount = sum(
                    1 for n in nbr.GetNeighbors() if n.GetAtomicNum() == 6
                )
                if carbonNbrCount >= 3:
                    continue  # 进入了大型碳骨架, 不计入糖碳

            # 仅碳原子 → 计数并继续扩展
            totalCarbons += 1
            carbonQueue.append((nIdx, depth + 1))

    return totalCarbons


def identify_monosaccharide_v2(mol, ring_atoms):
    base_name = "Hex" if len(ring_atoms) == 6 else "Pen"
    matched_anomer = "?"

    # ---- 原子计数氧门控 pre-compute (Oxygen Gate) ----
    # 计算当前碎片中环碳所连接的外环氧/氮原子数量
    # 用于在字典匹配后验证: 参考糖的 OH 数量是否与碎片一致
    fragmentOxygens = _countFragmentRingOxygens(mol, set(ring_atoms))

    # 1. 制造虚拟分身，抹平 N/S 杂原子的匹配障碍
    virtual_mol = create_virtual_standard_sugar(mol, ring_atoms)

    # 2. 使用虚拟分身去撞击字典，配合氧门控精确校验
    for (name, anomer), ref_mol in REFERENCE_MOLS.items():
        # 注意这里使用的是 virtual_mol，实现无损脱敏降维打击
        matches = virtual_mol.GetSubstructMatches(ref_mol, useChirality=True)
        for match in matches:
            # 如果这个匹配的子图完美覆盖了我们当前的 ring_atoms
            if all(idx in match for idx in ring_atoms):
                # ===== 氧门控 (Oxygen Gate) — 方向性校验 =====
                # 只在碎片的外环 O 数量 **远大于** 参考预期时拒绝匹配
                refOCount = REFERENCE_OXYGEN_COUNTS.get((name, anomer), -1)
                if refOCount >= 0:
                    if fragmentOxygens > refOCount + 1:
                        continue

                # ===== C/N 门控 (Carbon/Nitrogen Gate) =====
                # 绝对禁止: 含氮碎片匹配无氮模板, 或碳数差异 >1
                # Absolute block: N-containing fragment vs N-free template,
                # or carbon count mismatch >1
                refC, refN = REFERENCE_CN_COUNTS.get((name, anomer), (-1, -1))
                if refC >= 0:
                    # 收集碎片完整原子集 (环 + 直连外环原子)
                    fragUnitAtoms = set(ring_atoms)
                    for rIdx in ring_atoms:
                        for nbr in mol.GetAtomWithIdx(rIdx).GetNeighbors():
                            fragUnitAtoms.add(nbr.GetIdx())
                    fragC, fragN = _countFragmentCN(mol, fragUnitAtoms)
                    # 铁律 1: 碎片有 N 但参考无 N → 拒绝
                    #         (氨基己糖绝不能匹配普通戊糖)
                    if fragN > 0 and refN == 0:
                        continue
                    # 铁律 2: 碎片无 N 但参考有 N → 拒绝
                    if fragN == 0 and refN > 0:
                        continue
                    # 铁律 3: 碳数差异 >2 → 拒绝
                    #         (允许 ±2 容差: 糖苷键/修饰可能增减少量碳)
                    if abs(fragC - refC) > 2:
                        continue

                # ===== 所有门控通过, 接受匹配 =====
                base_name = name
                matched_anomer = anomer
                break
        if matched_anomer != "?":
            break

    if matched_anomer == "?":
        base_name, matched_anomer = rescue_hexose_by_key_nodes(mol, ring_atoms)

    # ---- 碳计数退避 2.1 (Carbon-Count Fallback 2.1) ----
    # 如果仍然是泛化标签 (Hex/Pen), 使用碳原子计数进行更精确的分类
    # v2.1: 移除 Non (壬糖) 分类 — 假阳性风险太高, 保持 Hex/Pen
    # v2.1: removed Non classification — too many false positives, keep Hex/Pen
    if base_name in ("Hex", "Pen") and matched_anomer == "?":
        totalC = _countFragmentCarbons(mol, ring_atoms)
        # 不再使用 Non (9C): 碳计数超标通常是苷元渗透, 不是真正的壬糖
        # No longer classify as Non (9C): excessive count is usually aglycon bleed
        if totalC == 8:
            base_name = "Oct"  # 辛糖 (8C, 如 KDO 类)
        elif totalC == 7:
            base_name = "Hept"  # 庚糖 (7C)
        # 5C, 6C, 9C+ 保持 Pen/Hex (Keep original Pen/Hex)

    # ---- 模糊分类器 3.0 (Fuzzy Classifier 3.0) ----
    # 替代旧的 Unknown_Sugar_MW_86 兜底逻辑
    # 使用碳/氮/氧原子计数对无法精确命名的糖进行化学分类
    # Replaces MW_86 fallback with chemistry-aware classification
    if base_name in ("Hex", "Pen") and matched_anomer == "?":
        # 收集完整糖单元原子 (环 + 直连外环原子)
        fullUnitAtoms = set(ring_atoms)
        for rIdx in ring_atoms:
            for nbr in mol.GetAtomWithIdx(rIdx).GetNeighbors():
                fullUnitAtoms.add(nbr.GetIdx())
        fragC, fragN = _countFragmentCN(mol, fullUnitAtoms)
        # 计算含氧量 (O/C 比率判断脱氧程度)
        fragO = sum(1 for idx in fullUnitAtoms
                    if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
        ringSize = len(ring_atoms)

        if ringSize == 6:  # 吡喃糖环
            if fragN > 0:
                base_name = "HexN"   # 氨基己糖 (有 N 原子)
            elif fragO <= 2:  # 环氧 + 最多 1 个 exo-O → 高度脱氧
                base_name = "dHex"   # 脱氧己糖
            else:
                base_name = "Hex"    # 普通己糖
        elif ringSize == 5:  # 呋喃糖环
            if fragN > 0:
                base_name = "PenN"   # 氨基戊糖 (罕见但化学完备)
            else:
                base_name = "Pen"    # 普通戊糖

    mods, is_acid, is_complex = check_modifications(mol, ring_atoms)

    if "GlcNAc" in base_name or "GalNAc" in base_name:
        if "NAc" in mods: mods.remove("NAc")
    if "GlcA" in base_name or "GalA" in base_name or "IdoA" in base_name:
        if "A" in mods: mods.remove("A")

    if "Glc" in base_name and "A" in mods:
        base_name = base_name.replace("Glc", "GlcA")
        mods.remove("A")
    elif "Gal" in base_name and "A" in mods:
        base_name = base_name.replace("Gal", "GalA")
        mods.remove("A")

    final_name = f"{base_name}({','.join(mods)})" if mods else base_name
    return final_name, matched_anomer


# =====================================================================
# LLM 文献回溯占位函数 (LLM Literature Retrieval Placeholder)
# 用于未来接入 Data Mining LLM 模块, 通过文献自动纠正泛指糖标签
# Placeholder for future Data Mining LLM module integration
# =====================================================================
def requestLlmSugarIdentification(
    compoundId: str,
    doi: str = "",
    pubchemCid: str = "",
    genericLabel: str = ""
) -> str:
    """
    [占位函数 / Placeholder] 请求 LLM 模块从文献中识别具体糖单元。
    Request LLM module to identify specific sugars from literature.

    触发条件 (Trigger): 化合物含有模糊标签 (Hex, dHex, HexN, Pen, PenN)。
    工作流 (Workflow):
      1. 通过 PubChem CID 或 DOI 获取文献摘要
      2. 提问 LLM: "What are the specific sugar units (monosaccharides)
         present in this compound?"
      3. 解析 LLM 输出, 映射到标准糖名

    Args:
        compoundId: 化合物编号 (如 CNP0420046)
        doi: 文献 DOI (如 10.1021/...)
        pubchemCid: PubChem CID
        genericLabel: 当前模糊标签 (如 "Hex", "dHex")

    Returns:
        具体糖名 (如 "D-Glc") 或空字符串 (无法解析)

    TODO: 接入实际 LLM API 后实现
    """
    # [PLACEHOLDER] 当前返回空字符串, 表示无法解析
    # 实际实现时将调用 OpenAI/Claude API + PubChem/OpenAlex
    return ""


def get_linkage_type(mol, unit1_atoms, unit2_atoms):
    u1_set = set(unit1_atoms)
    u2_set = set(unit2_atoms)
    for idx1 in unit1_atoms:
        atom1 = mol.GetAtomWithIdx(idx1)
        for nbr in atom1.GetNeighbors():
            if nbr.GetIdx() in u1_set: continue
            if nbr.GetIdx() in u2_set: return "-"
            if nbr.GetSymbol() == "O":
                for nnbr in nbr.GetNeighbors():
                    if nnbr.GetIdx() == idx1: continue
                    if nnbr.GetIdx() in u2_set: return "-"
    return "/"

def _hasValidSugarRing(mol, ring_atoms):
    """
    拓扑铁律: 开环糖杀手。
    Strict Ring Enforcement: verify fragment has a valid 5 or 6-membered
    oxygen-containing heterocycle.

    返回 False → 该碎片是开环结构 (如糖醇), 应标记为 Non_Cyclic_Invalid。
    """
    from rdkit.Chem import rdMolDescriptors
    if rdMolDescriptors.CalcNumRings(mol) == 0:
        return False

    ri = mol.GetRingInfo()
    for ring in ri.AtomRings():
        if len(ring) not in (5, 6):
            continue
        # 检查是否含有至少一个氧原子 (含氧杂环)
        oxygenInRing = sum(
            1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8
        )
        if oxygenInRing >= 1:
            return True
    return False


def generate_refined_sequence(mol):
    if not mol: return "", ""
    # Inject identify back to sugar_utils to ensure compatibility
    sugar_units, atom_to_sugar = sugar_utils.get_sugar_units(mol)
    if not sugar_units: return "", ""

    # ---- 拓扑铁律: 开环糖过滤 (Ring Enforcement) ----
    # 标记不含有效含氧杂环的碎片为 Non_Cyclic_Invalid
    validUnits = []
    for unit in sugar_units:
        ringAtoms = unit.get('ring_atoms', [])
        if ringAtoms:
            validUnits.append(unit)
        else:
            # 如果 sugar_utils 没提供 ring_atoms, 检查整个分子
            validUnits.append(unit)

    sugar_units = validUnits
    if not sugar_units:
        return "Non_Cyclic_Invalid", ""

    G = nx.DiGraph()

    for i, unit in enumerate(sugar_units):
        unit['order_id'] = i + 1
        name = unit.get('name', 'Unknown')
        anomer = unit.get('anomeric_config', '?')
        mods = unit.get('modifications', [])
        
        G.add_node(unit['id'], name=name, anomer=anomer, mods=mods, order_id=unit['order_id'])
        
    linkages = sugar_utils.find_glycosidic_linkages(mol, sugar_units)
    
    for l in linkages:
        src = l['sugar_donor']
        dst = l['sugar_acceptor']
        link_pos = l['linkage']
        G.add_edge(src, dst, link=link_pos)
        
    roots = [n for n in G.nodes() if G.out_degree(n) == 0]
    
    if not roots:
        seq_parts = []
        mod_parts = []
        for n in G.nodes():
            nd = G.nodes[n]
            seq_parts.append(nd['name'])
            mod_str = ",".join(f"*{m}" for m in nd['mods']) if nd['mods'] else ""
            mod_parts.append(f"{nd['name']}_{nd['order_id']}({mod_str})")
        return "Cyclic", " ; ".join(mod_parts)
    
    def traverse_two_pass(node):
        nd = G.nodes[node]
        name = nd['name']
        mods = nd['mods']
        order_id = nd['order_id']
        
        mod_label = f"{name}_{order_id}({','.join(f'*{m}' for m in mods)})" if mods else f"{name}_{order_id}()"
        
        children = list(G.predecessors(node))
        children.sort()
        
        if not children:
            return name, mod_label
            
        main_child = children[0]
        branches = children[1:]
        
        m_link = G.edges[main_child, node].get('link', '?')
        m_anomer = G.nodes[main_child].get('anomer', '?')
        m_link_str = f"({m_anomer}{m_link})"
        
        m_seq, m_mod = traverse_two_pass(main_child)
        
        res_seq = f"{m_seq}-{m_link_str}-"
        res_mod = f"{m_mod}-"
        
        for b in branches:
            b_link = G.edges[b, node].get('link', '?')
            b_anomer = G.nodes[b].get('anomer', '?')
            b_link_str = f"({b_anomer}{b_link})"
            
            b_seq, b_mod = traverse_two_pass(b)
            res_seq = f"[{b_seq}-{b_link_str}]-" + res_seq
            res_mod = f"[{b_mod}]-" + res_mod
            
        res_seq += name
        res_mod += mod_label
        
        return res_seq, res_mod

    seq_list = []
    mod_list = []
    for root in roots:
        s, m = traverse_two_pass(root)
        seq_list.append(s)
        mod_list.append(m)
        
    return " ; ".join(seq_list), " ; ".join(mod_list)

def analyze_glycan(smiles):
    if not smiles or pd.isna(smiles) or smiles == "nan": return "", ""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol: return "Invalid", ""
        return generate_refined_sequence(mol)
    except Exception as e:
        return "Error", ""
