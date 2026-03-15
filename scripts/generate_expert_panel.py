"""
GlycoNP Expert Review Panel — 500 样本可视化审查报告
GlycoNP Expert Review Panel — 500-Sample Visual Debug Report

为导师准备的专家审查面板: 每个化合物一张卡片,
左侧 2D 结构图 (四色高亮), 右侧化学/分类学参数。

三色+高亮 (Three-color + nucleotide highlight):
  红 = 糖核心骨架 (Sugar core)
  黄 = 修饰基团 (Modifications: O-Ac, NAc, etc.)
  绿 = 核苷酸碱基 (Nucleotide bases only — peptide detection REMOVED)
  蓝 = 苷元 (Aglycon, including lipid chains)

使用方法 (Usage):
  python scripts/generate_expert_panel.py [--input PATH] [--n 500]
"""
import argparse
import base64
import io
import os
import sys
import time
from typing import Dict, List, Optional, Set, Tuple

import pandas as pd
import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "lib")))


# =====================================================================
# 1. 采样逻辑 (Stratified Sampling)
# =====================================================================

def stratifiedSample(
    df: pd.DataFrame,
    nTotal: int = 500,
    nGlycolipids: int = 50,
    minClasses: int = 10,
) -> pd.DataFrame:
    """
    分层采样: 覆盖 10+ Superclasses, 保证 50 个 Glycolipids。
    Stratified sampling: cover 10+ Superclasses, guarantee 50 Glycolipids.

    策略 (Strategy):
      1. 先抽 nGlycolipids 个 Glycolipid
      2. 从剩余中按 Superclass 分配配额
      3. 确保至少 10 个不同大类被覆盖
    """
    df = df.copy()
    df["_class"] = df["Superclass"].apply(lambda x: _cleanClass(str(x)))

    # Step 1: 抽取 Glycolipids (Sample Glycolipids)
    glycolipidDf = df[df["_class"] == "Glycolipids"]
    nGL = min(nGlycolipids, len(glycolipidDf))
    sampledGL = glycolipidDf.sample(n=nGL, random_state=42) if nGL > 0 else pd.DataFrame()

    # Step 2: 稀有类保底 (Rare class minimum representation)
    remaining = df[~df.index.isin(sampledGL.index)]
    classCounts = remaining["_class"].value_counts()
    allClasses = classCounts.index.tolist()

    rareClasses = [c for c in allClasses if c != "Glycolipids"
                   and classCounts[c] < 20]

    rareSamples = []
    for cls in rareClasses[:max(0, minClasses - 2)]:
        clsDf = remaining[remaining["_class"] == cls]
        n = min(5, len(clsDf))
        if n > 0:
            rareSamples.append(clsDf.sample(n=n, random_state=42))

    rareDf = pd.concat(rareSamples) if rareSamples else pd.DataFrame()
    used = set(sampledGL.index) | set(rareDf.index)

    # Step 3: 剩余配额按比例分配 (Proportional fill)
    nRemaining = nTotal - len(sampledGL) - len(rareDf)
    pool = remaining[~remaining.index.isin(used)]

    if nRemaining > 0 and len(pool) > 0:
        fillDf = pool.sample(n=min(nRemaining, len(pool)), random_state=42)
    else:
        fillDf = pd.DataFrame()

    result = pd.concat([sampledGL, rareDf, fillDf], ignore_index=True)

    # 统计覆盖 (Coverage stats)
    coveredClasses = result["_class"].nunique()
    print(f"  Sampled: {len(result)} compounds")
    print(f"  Glycolipids: {len(result[result['_class'] == 'Glycolipids'])}")
    print(f"  Superclass coverage: {coveredClasses}")
    for cls in result["_class"].value_counts().head(12).items():
        print(f"    {cls[0]:30s} {cls[1]}")

    return result


def _cleanClass(val: str) -> str:
    if not val or val == "nan":
        return "Unclassified"
    s = str(val).strip()
    if "(Tanimoto=" in s:
        s = s[:s.index("(Tanimoto=")].strip()
    if s.startswith("Glycolipid"):
        s = "Glycolipids"
    return s if s else "Unclassified"


# =====================================================================
# 2. BFS 脂肪链修正 (Post-BFS Lipid Chain Correction)
# =====================================================================

def _reclassifyLipidChains(
    mol,
    glycanAtoms: Set[int],
    aglyconAtoms: Set[int],
    sugarUnits: List[Dict],
    minChainLength: int = 4,
) -> Tuple[Set[int], Set[int]]:
    """
    修正 BFS 洪水填充的误判: 脂肪链不应属于糖区。
    Fix BFS flood-fill misassignment: lipid chains should NOT be in glycan zone.

    问题原因 (Root Cause):
      BFS 从糖环出发, 通过酯键 (ester bond) 或醚键 (ether bond)
      漫入脂肪酸长链 → 整条碳链被误标为糖区 (红色)。

    修复策略 (Fix Strategy):
      1. 在糖区中找所有非环 sp3 碳原子
      2. 从每个这样的碳出发, DFS 沿非环碳链行走
      3. 如果连续非环碳链长度 > minChainLength, 整条链归苷元

    技术细节: 糖环上 C6-OH (如 CH2OH) 通常只有 1-2 个非环碳,
    不会被误切。只有长脂肪链 (>4 碳) 才会被重新归类。
    """
    from rdkit import Chem

    # 收集糖环原子 (Collect actual sugar ring atom indices)
    ringAtomSet: Set[int] = set()
    for unit in sugarUnits:
        for idx in unit.get("ring_atoms", []):
            ringAtomSet.add(idx)

    # 在糖区中找非环碳原子 (Find non-ring carbons in glycan zone)
    glycanNonRingCarbons = set()
    for idx in glycanAtoms:
        if idx in ringAtomSet:
            continue
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetAtomicNum() == 6 and not atom.IsInRing():
            glycanNonRingCarbons.add(idx)

    # DFS 找连续非环碳链 (DFS to find contiguous non-ring carbon chains)
    visited: Set[int] = set()
    chainsToReclassify: List[Set[int]] = []

    for startIdx in glycanNonRingCarbons:
        if startIdx in visited:
            continue

        # DFS 收集整条非环碳链 + 直接连接的杂原子
        chain: Set[int] = set()
        stack = [startIdx]
        while stack:
            cur = stack.pop()
            if cur in chain:
                continue
            atom = mol.GetAtomWithIdx(cur)
            # 只沿非环原子行走, 不进入糖环
            if cur in ringAtomSet:
                continue
            if atom.IsInRing():
                continue
            chain.add(cur)
            visited.add(cur)

            for nbr in atom.GetNeighbors():
                nbrIdx = nbr.GetIdx()
                if nbrIdx in chain or nbrIdx in ringAtomSet:
                    continue
                # 沿碳链行走 (Walk along carbon chain)
                if nbr.GetAtomicNum() == 6 and not nbr.IsInRing():
                    stack.append(nbrIdx)
                # 也收集链上直连的 O/N/S (如酯键O, 羟基O)
                elif nbr.GetAtomicNum() in (7, 8, 16) and not nbr.IsInRing():
                    # 但只收集度为 1-2 的杂原子 (不收集桥O)
                    if nbr.GetDegree() <= 2:
                        chain.add(nbrIdx)

        # 计算链中碳原子数 (Count carbons in chain)
        carbonCount = sum(
            1 for idx in chain
            if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6
        )

        if carbonCount > minChainLength:
            chainsToReclassify.append(chain)

    # 重新归类 (Reclassify)
    totalReclassified = 0
    for chain in chainsToReclassify:
        for idx in chain:
            if idx in glycanAtoms:
                glycanAtoms.discard(idx)
                aglyconAtoms.add(idx)
                totalReclassified += 1

    return glycanAtoms, aglyconAtoms


# =====================================================================
# 3. 分子图生成 (Molecule Image Rendering)
# =====================================================================

def renderMoleculeImage(
    smiles: str,
    imgSize: Tuple[int, int] = (500, 350),
) -> Optional[str]:
    """
    生成四色高亮的分子 PNG, 返回 base64 编码。
    Generate four-color highlighted molecule PNG, return as base64.

    四色 (Four colors):
      红 (1.0, 0.70, 0.70) = 糖核心
      黄 (1.0, 0.90, 0.50) = 修饰基团
      绿 (0.70, 1.0, 0.80)  = 核苷酸/肽键
      蓝 (0.70, 0.80, 1.0)  = 苷元
    """
    from rdkit import Chem
    from rdkit.Chem import Draw, AllChem

    if not smiles or str(smiles) in ("nan", "", "NULL", "*"):
        return None

    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        return None

    AllChem.Compute2DCoords(mol)

    # 尝试划分糖区 (Try zone assignment)
    try:
        from sugar_utils import find_mapped_sugar_units
        sugarUnits = find_mapped_sugar_units(mol)
    except Exception:
        sugarUnits = []

    if not sugarUnits:
        # 无糖单元: 全部蓝色或无高亮 (No sugar: plain drawing)
        drawer = Draw.rdMolDraw2D.MolDraw2DCairo(imgSize[0], imgSize[1])
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        pngData = drawer.GetDrawingText()
        return base64.b64encode(pngData).decode("ascii")

    # 区域划分 (Zone assignment)
    try:
        from phase7_visualizer import identifySugarAtomZones
        glycanAtoms, aglyconAtoms = identifySugarAtomZones(mol, sugarUnits)
    except Exception:
        glycanAtoms = set()
        aglyconAtoms = set(range(mol.GetNumAtoms()))

    # BFS 修正: 将误入糖区的脂肪链重新归类为苷元
    # Post-BFS correction: reclassify lipid chains that leaked into glycan zone
    glycanAtoms, aglyconAtoms = _reclassifyLipidChains(
        mol, glycanAtoms, aglyconAtoms, sugarUnits)

    # 修饰原子检测 (Modification atoms)
    modAtoms: Set[int] = set()
    try:
        from glycan_modifications import _getCompiledModSmarts
        smartsDict = _getCompiledModSmarts()
        for modName, pattern in smartsDict.items():
            if pattern is None:
                continue
            matches = mol.GetSubstructMatches(pattern)
            for match in matches:
                matchSet = set(match)
                if matchSet & glycanAtoms:
                    modAtoms.update(matchSet & glycanAtoms)
    except Exception:
        pass

    # 核苷酸碱基检测 — 仅保留嘌呤/嘧啶 (Nucleotide bases only)
    # 已移除肽键 SMARTS: [NX3][CX3](=[OX1])[CX4] 会误判 NAc 为氨基酸
    # REMOVED peptide SMARTS: it falsely matches NAc (N-Acetyl) on sugars
    greenAtoms: Set[int] = set()
    try:
        # 嘌呤碱基 (Purine base — adenine/guanine)
        purine = Chem.MolFromSmarts("[#7]1~[#6]~[#7]~[#6]2~[#6]1~[#7]~[#6]~[#7]~2")
        # 嘧啶碱基 (Pyrimidine base — cytosine/uracil/thymine)
        pyrimidine = Chem.MolFromSmarts("[#7]1~[#6]~[#6]~[#7]~[#6](~[#8])~[#7]~1")

        for pat in [purine, pyrimidine]:
            if pat is not None:
                for match in mol.GetSubstructMatches(pat):
                    greenAtoms.update(match)
    except Exception:
        pass

    # 四色分配 (Four-color assignment)
    COLOR_SUGAR = (1.0, 0.70, 0.70)
    COLOR_MOD = (1.0, 0.90, 0.50)
    COLOR_GREEN = (0.70, 1.0, 0.80)
    COLOR_AGLYCON = (0.70, 0.80, 1.0)

    atomColorMap = {}
    bondColorMap = {}

    for idx in range(mol.GetNumAtoms()):
        if idx in greenAtoms:
            atomColorMap[idx] = COLOR_GREEN
        elif idx in modAtoms:
            atomColorMap[idx] = COLOR_MOD
        elif idx in glycanAtoms:
            atomColorMap[idx] = COLOR_SUGAR
        elif idx in aglyconAtoms:
            atomColorMap[idx] = COLOR_AGLYCON

    for bond in mol.GetBonds():
        a = bond.GetBeginAtomIdx()
        b = bond.GetEndAtomIdx()
        if a in greenAtoms and b in greenAtoms:
            bondColorMap[bond.GetIdx()] = COLOR_GREEN
        elif a in modAtoms and b in modAtoms:
            bondColorMap[bond.GetIdx()] = COLOR_MOD
        elif a in glycanAtoms and b in glycanAtoms:
            bondColorMap[bond.GetIdx()] = COLOR_SUGAR
        elif a in aglyconAtoms and b in aglyconAtoms:
            bondColorMap[bond.GetIdx()] = COLOR_AGLYCON

    drawer = Draw.rdMolDraw2D.MolDraw2DCairo(imgSize[0], imgSize[1])
    drawer.DrawMolecule(
        mol,
        highlightAtoms=list(atomColorMap.keys()),
        highlightAtomColors=atomColorMap,
        highlightBonds=list(bondColorMap.keys()),
        highlightBondColors=bondColorMap,
    )
    drawer.FinishDrawing()
    pngData = drawer.GetDrawingText()
    return base64.b64encode(pngData).decode("ascii")


# =====================================================================
# 3. HTML 报告生成 (HTML Report Generation)
# =====================================================================

def generateHtmlReport(
    sampleDf: pd.DataFrame,
    outputPath: str,
):
    """
    生成专家审查面板 HTML: 左图右参。
    Generate Expert Review Panel HTML: left=image, right=parameters.
    """

    htmlParts = []
    htmlParts.append(_htmlHeader())

    total = len(sampleDf)
    success = 0
    failed = 0

    for i, (_, row) in enumerate(sampleDf.iterrows()):
        if i % 50 == 0:
            print(f"  Rendering {i+1}/{total}...", end="\r")

        smiles = str(row.get("canonical_smiles", ""))
        imgB64 = renderMoleculeImage(smiles)
        if imgB64:
            success += 1
        else:
            failed += 1

        htmlParts.append(_htmlCard(row, imgB64, i + 1))

    print(f"\n  Rendered: {success} success, {failed} failed")

    htmlParts.append(_htmlFooter())

    with open(outputPath, "w", encoding="utf-8") as f:
        f.write("\n".join(htmlParts))

    print(f"  HTML saved: {outputPath}")


def _htmlHeader() -> str:
    return """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>GlycoNP Expert Review Panel</title>
<style>
  * { margin: 0; padding: 0; box-sizing: border-box; }
  body {
    font-family: 'Segoe UI', 'Helvetica Neue', Arial, sans-serif;
    background: #0f0f1a;
    color: #e0e0e0;
    padding: 20px;
  }
  h1 {
    text-align: center;
    color: #fff;
    margin-bottom: 10px;
    font-size: 28px;
  }
  .subtitle {
    text-align: center;
    color: #888;
    margin-bottom: 30px;
    font-size: 14px;
  }
  .legend {
    display: flex;
    justify-content: center;
    gap: 20px;
    margin-bottom: 30px;
    flex-wrap: wrap;
  }
  .legend-item {
    display: flex;
    align-items: center;
    gap: 6px;
    font-size: 13px;
  }
  .legend-dot {
    width: 14px; height: 14px;
    border-radius: 50%;
    display: inline-block;
  }
  .card {
    display: flex;
    background: #1a1a2e;
    border: 1px solid #2a2a4a;
    border-radius: 12px;
    margin-bottom: 16px;
    overflow: hidden;
    transition: border-color 0.2s;
  }
  .card:hover {
    border-color: #4a90d9;
  }
  .card-img {
    flex: 0 0 520px;
    padding: 12px;
    background: #fff;
    display: flex;
    flex-direction: column;
    align-items: center;
    justify-content: center;
  }
  .card-img img {
    max-width: 500px;
    max-height: 350px;
  }
  .card-img .no-img {
    color: #999;
    font-size: 14px;
  }
  .card-info {
    flex: 1;
    padding: 16px 20px;
    display: flex;
    flex-direction: column;
    gap: 6px;
    font-size: 13px;
    line-height: 1.5;
  }
  .card-id {
    font-size: 18px;
    font-weight: 700;
    color: #4a90d9;
    margin-bottom: 4px;
  }
  .card-name {
    font-size: 15px;
    color: #f5a623;
    margin-bottom: 8px;
    word-break: break-word;
  }
  .info-row {
    display: flex;
    gap: 8px;
  }
  .info-label {
    color: #888;
    min-width: 130px;
    flex-shrink: 0;
  }
  .info-value {
    color: #ddd;
    word-break: break-all;
  }
  .info-value.sugar {
    color: #ff8888;
    font-weight: 600;
  }
  .info-value.mod {
    color: #ffdd66;
  }
  .info-value.class-tag {
    color: #88ccff;
    font-weight: 600;
  }
  .badge {
    display: inline-block;
    padding: 2px 8px;
    border-radius: 10px;
    font-size: 11px;
    font-weight: 600;
    margin-right: 4px;
  }
  .badge-nuc { background: #2d5a3d; color: #7fdf9f; }
  .badge-pep { background: #5a3d2d; color: #df9f7f; }
  .badge-gl  { background: #2d3d5a; color: #7f9fdf; }
  a { color: #4a90d9; text-decoration: none; }
  a:hover { text-decoration: underline; }
  .top-banner {
    text-align: center;
    font-size: 12px;
    color: #666;
    padding: 8px;
    margin-bottom: 6px;
  }
  .bot-banner {
    text-align: center;
    font-size: 12px;
    color: #aaa;
    padding: 4px;
    background: rgba(255,255,255,0.05);
    border-radius: 6px;
    margin-top: 4px;
  }
</style>
</head>
<body>
<h1>GlycoNP Expert Review Panel</h1>
<p class="subtitle">500-Sample Visual Debug Report | Pipeline v7</p>
<div class="legend">
  <div class="legend-item">
    <span class="legend-dot" style="background:#ff8888"></span> Sugar Core
  </div>
  <div class="legend-item">
    <span class="legend-dot" style="background:#ffdd66"></span> Modifications (O-Ac, NAc, etc.)
  </div>
  <div class="legend-item">
    <span class="legend-dot" style="background:#88ffcc"></span> Nucleotide Base
  </div>
  <div class="legend-item">
    <span class="legend-dot" style="background:#88aaff"></span> Aglycon (incl. lipid chains)
  </div>
</div>
"""


def _htmlCard(row: pd.Series, imgB64: Optional[str], idx: int) -> str:
    """生成单个化合物卡片 HTML."""
    name = str(row.get("name", ""))[:60]
    if name == "nan":
        name = ""
    inchiKey = str(row.get("standard_inchi_key", ""))
    superclass = _cleanClass(str(row.get("Superclass", "")))
    sugar = str(row.get("Sugar_Sequence", ""))
    if sugar == "nan":
        sugar = ""
    mods = str(row.get("Glycan_Modifications", ""))
    if mods == "nan":
        mods = "None"
    organism = str(row.get("organisms", ""))
    if organism == "nan":
        organism = "Unknown"
    else:
        organism = organism.split("|")[0][:50]
    scaffold = str(row.get("Murcko_Scaffold", ""))[:50]
    if scaffold == "nan":
        scaffold = ""
    hasNuc = str(row.get("Has_Nucleotide", "")) == "True"
    hasPep = str(row.get("Has_Peptide", "")) == "True"
    dois = str(row.get("dois", ""))

    # 分子量等
    mw = str(row.get("molecular_weight", ""))
    alogp = str(row.get("alogp", ""))
    npLike = str(row.get("np_likeness", ""))
    formula = str(row.get("molecular_formula", ""))

    # DOI 链接 (DOI links)
    doiHtml = ""
    if dois and dois != "nan":
        doiList = dois.split("|")[:3]
        doiLinks = []
        for d in doiList:
            d = d.strip()
            if d:
                url = d if d.startswith("http") else f"https://doi.org/{d}"
                doiLinks.append(f'<a href="{url}" target="_blank">{d[:40]}</a>')
        if doiLinks:
            doiHtml = " | ".join(doiLinks)

    # Badges
    badges = ""
    if hasNuc:
        badges += '<span class="badge badge-nuc">Nucleotide</span>'
    if hasPep:
        badges += '<span class="badge badge-pep">Peptide</span>'
    if superclass == "Glycolipids":
        badges += '<span class="badge badge-gl">Glycolipid</span>'

    # 图片 (Image)
    if imgB64:
        imgHtml = f'<img src="data:image/png;base64,{imgB64}" alt="Structure">'
    else:
        imgHtml = '<span class="no-img">Structure unavailable</span>'

    # 图上下方标注 (Top/bottom annotations on image)
    topBanner = f"{organism} | {superclass}"
    botBanner = f"{sugar}" + (f" | {mods}" if mods != "None" else "")

    return f"""
<div class="card">
  <div class="card-img">
    <div class="top-banner">{topBanner}</div>
    {imgHtml}
    <div class="bot-banner">{botBanner}</div>
  </div>
  <div class="card-info">
    <div class="card-id">#{idx} {badges}</div>
    <div class="card-name">{name}</div>
    <div class="info-row"><span class="info-label">InChIKey</span><span class="info-value">{inchiKey}</span></div>
    <div class="info-row"><span class="info-label">Superclass</span><span class="info-value class-tag">{superclass}</span></div>
    <div class="info-row"><span class="info-label">Sugar Sequence</span><span class="info-value sugar">{sugar}</span></div>
    <div class="info-row"><span class="info-label">Modifications</span><span class="info-value mod">{mods}</span></div>
    <div class="info-row"><span class="info-label">Organism</span><span class="info-value">{organism}</span></div>
    <div class="info-row"><span class="info-label">Murcko Scaffold</span><span class="info-value" style="font-size:11px">{scaffold}</span></div>
    <div class="info-row"><span class="info-label">Mol Formula</span><span class="info-value">{formula}</span></div>
    <div class="info-row"><span class="info-label">MW / ALogP</span><span class="info-value">{mw} / {alogp}</span></div>
    <div class="info-row"><span class="info-label">NP-Likeness</span><span class="info-value">{npLike}</span></div>
    <div class="info-row"><span class="info-label">Has Nucleotide</span><span class="info-value">{"Yes" if hasNuc else "No"}</span></div>
    <div class="info-row"><span class="info-label">Has Peptide</span><span class="info-value">{"Yes" if hasPep else "No"}</span></div>
    <div class="info-row"><span class="info-label">Literature</span><span class="info-value">{doiHtml if doiHtml else "No DOI available"}</span></div>
  </div>
</div>"""


def _htmlFooter() -> str:
    return """
<p style="text-align:center; color:#555; margin-top:40px; font-size:12px;">
  Generated by GlycoNP Pipeline | Expert Review Panel
</p>
</body>
</html>"""


# =====================================================================
# Main
# =====================================================================

def main():
    parser = argparse.ArgumentParser(
        description="GlycoNP Expert Review Panel")
    parser.add_argument("--input", type=str, default=None, help="Input CSV")
    parser.add_argument("--n", type=int, default=500, help="Number of samples")
    parser.add_argument("--glycolipids", type=int, default=50,
                        help="Number of guaranteed Glycolipids")
    args = parser.parse_args()

    baseDir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    reportDir = os.path.join(baseDir, "reports")

    if args.input:
        inputPath = args.input
    else:
        inputPath = os.path.join(reportDir, "GlycoNP_Pipeline_Full.csv")
        if not os.path.exists(inputPath):
            inputPath = os.path.join(reportDir, "GlycoNP_Pipeline_1000.csv")

    print("=" * 60)
    print("  GlycoNP Expert Review Panel")
    print("=" * 60)
    print(f"  Input: {inputPath}")

    df = pd.read_csv(inputPath, low_memory=False, dtype=str, encoding="utf-8-sig")
    print(f"  Loaded {len(df):,} rows")

    sampleDf = stratifiedSample(df, nTotal=args.n, nGlycolipids=args.glycolipids)

    outputPath = os.path.join(reportDir, "Expert_Review_Panel.html")
    t0 = time.time()
    generateHtmlReport(sampleDf, outputPath)
    elapsed = time.time() - t0

    print(f"\n{'='*60}")
    print(f"  Done! ({elapsed:.0f}s)")
    print(f"  Output: {outputPath}")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
