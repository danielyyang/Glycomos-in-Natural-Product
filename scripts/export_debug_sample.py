"""
GlycoNP Colored Debug Sample Exporter — 500 Compounds with 3-Color Highlighting
(三色高亮调试样本导出器)

Full Molecule column uses phase7_visualizer's drawHighlightedMolecule():
  Red    = Sugar core skeleton (糖核心骨架)
  Yellow = Modification groups (修饰基团: O-Ac, NAc, Sulfate)
  Blue   = Aglycone / non-sugar (苷元)

Usage / 用法:
  python scripts/export_debug_sample.py
"""
import base64
import io
import os
import sys
import time

import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Draw, AllChem

# 确保 lib/ 目录在路径中 (Add lib/ to path)
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "lib")))

SAMPLE_SIZE = 500
IMG_SIZE = (350, 260)
IMG_SIZE_SMALL = (280, 200)

# =====================================================================
# 图像生成函数 (Image Generation Functions)
# =====================================================================

def coloredMolToBase64(smiles: str, size: tuple = IMG_SIZE) -> tuple:
    """
    使用 phase7_visualizer 的三色高亮渲染完整分子。
    Render full molecule with three-color highlighting via phase7_visualizer.

    红色 = 糖核心, 黄色 = 修饰基团, 蓝色 = 苷元
    Red = Sugar core, Yellow = Modifications, Blue = Aglycone

    Returns:
        (imgHtml: str, highlighted: bool)
    """
    if not smiles or str(smiles) in ("nan", "", "None", "Aliphatic Chain"):
        return '<span style="color:#aaa;font-size:11px;">N/A</span>', False
    try:
        from phase7_visualizer import drawHighlightedMolecule
        pngBytes = drawHighlightedMolecule(
            smiles=str(smiles),
            sugarUnits=None,  # 自动检测 (auto-detect)
            imgSize=size,
        )
        if pngBytes is None:
            return plainMolToBase64(smiles, size), False
        b64 = base64.b64encode(pngBytes).decode("utf-8")
        imgTag = f'<img src="data:image/png;base64,{b64}" width="{size[0]}" height="{size[1]}" />'
        return imgTag, True
    except Exception as e:
        # 回退到无高亮渲染 (Fallback to plain rendering)
        return plainMolToBase64(smiles, size), False


def plainMolToBase64(smiles: str, size: tuple = IMG_SIZE_SMALL) -> str:
    """
    普通无高亮 SMILES 渲染 (用于 Aglycon/Glycan 单独列)。
    Plain SMILES rendering without highlighting (for Aglycon/Glycan columns).
    Handles dummy atoms from MolFragmentToSmiles ([200*], [100*], etc.)
    """
    if not smiles or str(smiles) in ("nan", "", "None", "Aliphatic Chain"):
        return '<span style="color:#aaa;font-size:11px;">N/A</span>'
    try:
        import re
        smi = str(smiles)
        # 清理 MolFragmentToSmiles 产生的编号 dummy 原子
        # Strip numbered dummy atoms: [200*], [100*] → [*]
        smi = re.sub(r'\[\d+\*\]', '[*]', smi)
        # 尝试解析 (Try parsing with sanitization)
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            # 降级尝试: 不 sanitize (Fallback: no sanitization)
            mol = Chem.MolFromSmiles(smi, sanitize=False)
        if mol is None:
            return '<span style="color:#c00;font-size:11px;">Parse Error</span>'
        AllChem.Compute2DCoords(mol)
        img = Draw.MolToImage(mol, size=size)
        buf = io.BytesIO()
        img.save(buf, format="PNG")
        b64 = base64.b64encode(buf.getvalue()).decode("utf-8")
        return f'<img src="data:image/png;base64,{b64}" width="{size[0]}" height="{size[1]}" />'
    except Exception as e:
        return f'<span style="color:#c00;font-size:10px;">Error: {str(e)[:30]}</span>'


def main():
    baseDir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    reportDir = os.path.join(baseDir, "reports")
    inputPath = os.path.join(reportDir, "GlycoNP_Pipeline_Full_Cleaned.csv")

    print("=" * 70)
    print("  GlycoNP Colored Debug Sample — 500 Visual HTML")
    print("  三色高亮: Red=糖骨架 Yellow=修饰 Blue=苷元")
    print("=" * 70)

    t0 = time.time()
    df = pd.read_csv(inputPath, low_memory=False, encoding="utf-8-sig", dtype=str)
    total = len(df)
    print(f"  Loaded: {total:,} rows")

    # ---- MW_86 调查 (MW_86 Investigation) ----
    mw86Df = df[df["Sugar_Sequence"].str.contains("Unknown_Sugar_MW_86", na=False)]
    print(f"\n  === MW_86 Investigation ===")
    print(f"  Total Unknown_Sugar_MW_86 rows: {len(mw86Df):,}")
    print(f"  Sample 10 SMILES:")
    for i, (_, row) in enumerate(mw86Df.head(10).iterrows()):
        cid = str(row.get("identifier", "?"))[:20]
        glySmiles = str(row.get("Glycan_SMILES", ""))[:100]
        seq = str(row.get("Sugar_Sequence", ""))[:70]
        print(f"    {i+1}. [{cid}] Seq={seq}")
        print(f"       Glycan={glySmiles}")
    print()

    # ---- 过滤有效行 (Filter valid rows) ----
    validDf = df[
        df["canonical_smiles"].notna()
        & (df["canonical_smiles"] != "")
        & df["Sugar_Sequence"].notna()
        & (df["Sugar_Sequence"] != "")
        & (df["Sugar_Sequence"] != "nan")
    ].copy()
    print(f"  Valid rows: {len(validDf):,}")

    # 强制包含目标化合物 (Force-include target compounds for debug verification)
    TARGET_IDS = ["CNP0426561", "CNP0207797"]
    targetDf = validDf[validDf["identifier"].str.contains("|".join(TARGET_IDS), na=False)]
    print(f"  Force-included target compounds: {len(targetDf)} rows")

    # 随机抽样 (Random sample from remaining)
    np.random.seed(42)
    remainDf = validDf[~validDf.index.isin(targetDf.index)]
    randomDf = remainDf.sample(n=min(SAMPLE_SIZE - len(targetDf), len(remainDf)), random_state=42)
    sampleDf = pd.concat([targetDf, randomDf]).reset_index(drop=True)
    print(f"  Sampled: {len(sampleDf)} rows (including {len(targetDf)} targets)")

    # ---- 生成图片 (Generate images) ----
    print(f"  Generating RDKit 2D images with 3-color highlighting...")
    rows = []
    highlightSuccessCount = 0
    highlightFailCount = 0

    for i, (idx, row) in enumerate(sampleDf.iterrows()):
        if (i + 1) % 50 == 0:
            print(f"    {i+1}/{len(sampleDf)}...")

        compoundId = str(row.get("identifier", ""))[:20]
        compoundName = str(row.get("name", ""))[:60]
        fullSmiles = str(row.get("canonical_smiles", ""))
        sugarSeqRaw = str(row.get("Sugar_Sequence", ""))
        organism = str(row.get("organisms", ""))[:80]

        # === 实时 Sugar Sequence 重计算, 消除 CSV 中残留的 Non ===
        # Live recomputation to eliminate stale 'Non' from old CSV data
        sugarSeq = sugarSeqRaw
        try:
            liveMolSeq = Chem.MolFromSmiles(str(fullSmiles))
            if liveMolSeq:
                from sugar_utils import find_mapped_sugar_units
                liveUnitsSeq = find_mapped_sugar_units(liveMolSeq)
                if liveUnitsSeq:
                    sugarSeq = " → ".join(u["name"] for u in liveUnitsSeq)
        except Exception:
            pass  # 回退到 CSV 数据 (Fallback to CSV data)

        # === 实时修饰重计算 v2.0 (Live Modification Recomputation v2.0) ===
        # 使用新的直接扫描器重新检测修饰, 而非读取 CSV 中的旧数据
        # Re-detect modifications using the new v2.0 direct scanner
        mods = ""
        try:
            liveMol = Chem.MolFromSmiles(str(fullSmiles))
            if liveMol:
                from sugar_utils import find_mapped_sugar_units
                liveUnits = find_mapped_sugar_units(liveMol)
                modParts = []
                for u in liveUnits:
                    uMods = u.get('modifications', [])
                    if uMods:
                        modParts.append(f"{u['name']}({','.join('*'+m for m in uMods)})")
                    else:
                        modParts.append(f"{u['name']}()")
                mods = " ; ".join(modParts)
        except Exception:
            mods = str(row.get("Glycan_Modifications", ""))
        superclass = str(row.get("Superclass", ""))
        if "(Tanimoto=" in superclass:
            superclass = superclass[:superclass.index("(Tanimoto=")].strip()

        # 完整分子用三色高亮 (Full molecule with 3-color highlight)
        fullImg, wasHighlighted = coloredMolToBase64(fullSmiles, IMG_SIZE)
        if wasHighlighted:
            highlightSuccessCount += 1
        else:
            highlightFailCount += 1

        # ============================================================
        # 动态提取苷元/糖链片段 — 使用键切断+封端 (Bond Cleavage + Capping)
        # Uses FragmentOnBonds() for chemically valid fragment extraction
        # ============================================================
        aglyImg = '<span style="color:#aaa;font-size:11px;">N/A</span>'
        glyImg = '<span style="color:#aaa;font-size:11px;">N/A</span>'
        try:
            mol = Chem.MolFromSmiles(str(fullSmiles))
            if mol is not None:
                from sugar_utils import find_mapped_sugar_units
                from phase7_visualizer import identifySugarAtomZones, cleaveAndCap
                sugarUnitsLocal = find_mapped_sugar_units(mol)
                if sugarUnitsLocal:
                    # v2.1: identifySugarAtomZones 返回 3 个集合
                    # v2.1: returns 3 sets (ring, substituent, aglycon)
                    sugarRingAtoms, substituentAtoms, aglyconAtomSet = identifySugarAtomZones(mol, sugarUnitsLocal)
                    glycanAtomSet = sugarRingAtoms | substituentAtoms
                    # 使用键切断+封端 (Use bond cleavage + capping)
                    aglyconSmi, glycanSmi = cleaveAndCap(
                        mol, glycanAtomSet, aglyconAtomSet,
                        minAglyconHeavyAtoms=3,
                    )
                    if aglyconSmi:
                        aglyImg = plainMolToBase64(aglyconSmi, IMG_SIZE_SMALL)
                    else:
                        aglyImg = '<span style="color:#aaa;font-size:11px;">No Aglycone</span>'
                    if glycanSmi:
                        glyImg = plainMolToBase64(glycanSmi, IMG_SIZE_SMALL)
        except Exception:
            pass

        # === 实时 Name-Rescue 重计算 (Live Name-Rescue Recomputation) ===
        # 使用修复后的分号分割+位置映射逻辑重新计算
        rescuedSeq = ""
        rescueMethod = ""
        try:
            from stereo_rescue import rescueSugarSequence
            compName = str(row.get("name", ""))
            iupacName = str(row.get("iupac_name", ""))
            rescuedSeq, rescueMethod = rescueSugarSequence(
                smiles=fullSmiles,
                currentSequence=sugarSeq,
                name=compName if compName != "nan" else None,
                iupacName=iupacName if iupacName != "nan" else None,
            )
            if rescuedSeq == sugarSeq:
                rescuedSeq = ""
                rescueMethod = ""
        except Exception:
            rescuedSeq = str(row.get("Rescued_Sugar_Sequence", ""))
            rescueMethod = str(row.get("Rescue_Method", ""))

        rows.append({
            "ID": compoundId,
            "Name": compoundName if compoundName != "nan" else "",
            "Full Molecule (Colored)": fullImg,
            "Aglycon": aglyImg,
            "Glycan": glyImg,
            "Sugar Sequence": sugarSeq,
            "Rescued Seq": rescuedSeq if rescuedSeq not in ("nan", "") else "",
            "Rescue": rescueMethod if rescueMethod not in ("nan", "") else "",
            "Modifications": mods if mods != "nan" else "",
            "Organism": organism if organism != "nan" else "",
            "Superclass": superclass if superclass != "nan" else "",
        })

    htmlDf = pd.DataFrame(rows)

    # ---- 生成 HTML (Build HTML) ----
    print(f"  Building HTML with color legend...")

    css = """
    <style>
        body { font-family: 'Segoe UI', Arial, sans-serif; margin: 20px; background: #f8f9fa; }
        h1 { color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 10px; }
        .stats { background: #ecf0f1; padding: 12px 20px; border-radius: 8px; margin-bottom: 10px;
                 font-size: 14px; color: #34495e; }
        .legend { background: white; padding: 12px 20px; border-radius: 8px; margin-bottom: 20px;
                  font-size: 13px; border: 2px solid #3498db; display: flex; gap: 30px; align-items: center; }
        .legend-item { display: flex; align-items: center; gap: 6px; }
        .legend-swatch { width: 20px; height: 20px; border-radius: 4px; border: 1px solid #aaa; }
        table { border-collapse: collapse; width: 100%; background: white; box-shadow: 0 2px 8px rgba(0,0,0,0.1); }
        th { background: #2c3e50; color: white; padding: 10px 6px; font-size: 12px;
             position: sticky; top: 0; z-index: 10; }
        td { border: 1px solid #ddd; padding: 6px; font-size: 11px; vertical-align: middle;
             max-width: 220px; word-wrap: break-word; }
        tr:nth-child(even) { background: #f5f6fa; }
        tr:hover { background: #eaf2f8; }
        img { border: 1px solid #ddd; border-radius: 4px; }
        .seq { font-family: 'Consolas', monospace; font-size: 11px; color: #c0392b; font-weight: bold; }
        .org { font-size: 10px; color: #7f8c8d; }
        .cls { font-size: 11px; color: #2980b9; font-weight: bold; }
        .mod { font-size: 10px; color: #8e44ad; }
    </style>
    """

    # 颜色图例 (Color legend)
    legend = """
    <div class="legend">
        <b>🎨 Color Legend / 颜色图例:</b>
        <div class="legend-item">
            <div class="legend-swatch" style="background: rgb(255,179,179);"></div>
            <span>Sugar Core (糖核心骨架)</span>
        </div>
        <div class="legend-item">
            <div class="legend-swatch" style="background: rgb(255,230,128);"></div>
            <span>Modifications (修饰基团)</span>
        </div>
        <div class="legend-item">
            <div class="legend-swatch" style="background: rgb(179,204,255);"></div>
            <span>Aglycone (苷元)</span>
        </div>
    </div>
    """

    def formatCell(col, val):
        if col == "Sugar Sequence":
            return f'<span class="seq">{val}</span>'
        elif col == "Organism":
            return f'<span class="org">{val}</span>'
        elif col == "Superclass":
            return f'<span class="cls">{val}</span>'
        elif col == "Modifications":
            return f'<span class="mod">{val}</span>'
        return val

    tableHtml = '<table>\n<thead><tr>'
    for col in htmlDf.columns:
        tableHtml += f'<th>{col}</th>'
    tableHtml += '</tr></thead>\n<tbody>\n'

    for _, row in htmlDf.iterrows():
        tableHtml += '<tr>'
        for col in htmlDf.columns:
            val = row[col]
            if col in ("Full Molecule (Colored)", "Aglycon", "Glycan"):
                tableHtml += f'<td>{val}</td>'
            else:
                tableHtml += f'<td>{formatCell(col, val)}</td>'
        tableHtml += '</tr>\n'
    tableHtml += '</tbody></table>'

    fullHtml = f"""<!DOCTYPE html>
<html><head>
<meta charset="utf-8">
<title>GlycoNP Colored Debug Sample — {SAMPLE_SIZE} Compounds</title>
{css}
</head><body>
<h1>🔬 GlycoNP Colored Debug Sample — {len(sampleDf)} Compounds</h1>
<div class="stats">
  <b>Dataset:</b> GlycoNP_Pipeline_Full_Cleaned.csv |
  <b>Total:</b> {total:,} |
  <b>Sampled:</b> {len(sampleDf)} (random, seed=42) |
  <b>Generated:</b> {time.strftime('%Y-%m-%d %H:%M')} |
  <b>Highlight Success:</b> {highlightSuccessCount} |
  <b>Highlight Fallback:</b> {highlightFailCount}
</div>
{legend}
{tableHtml}
</body></html>"""

    outPath = os.path.join(reportDir, "debug_sample_500_colored.html")
    with open(outPath, "w", encoding="utf-8") as f:
        f.write(fullHtml)

    print(f"\n  Output: {outPath}")
    print(f"  Size: {os.path.getsize(outPath)/1024/1024:.1f} MB")
    print(f"  Highlight success: {highlightSuccessCount}")
    print(f"  Highlight fallback: {highlightFailCount}")
    print(f"  Time: {time.time()-t0:.0f}s")
    print(f"  {'='*70}")


if __name__ == "__main__":
    main()
