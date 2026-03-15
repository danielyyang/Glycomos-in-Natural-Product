"""
GlycoNP UMAP 化学空间聚类可视化
GlycoNP UMAP Chemical Space Clustering

将苷元 (Aglycon) 的 Morgan 指纹通过 UMAP 降维到 2D, 生成交互式散点图。
直接回答导师的终极问题:
  "在苷元宇宙中, 特定糖链是否像云团一样只精准地挂在特定的苷元岛屿上?"

使用方法 (Usage):
  python scripts/plot_umap_chemical_space.py [--input PATH] [--sample 10000]
"""
import argparse
import os
import sys
import time

import numpy as np
import pandas as pd

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))


# =====================================================================
# 1. Morgan 指纹计算 (Morgan Fingerprint Computation)
# =====================================================================

def computeMorganFingerprints(
    smilesList: list,
    radius: int = 2,
    nBits: int = 2048,
) -> np.ndarray:
    """
    从 SMILES 列表计算 Morgan 指纹矩阵。
    Compute Morgan fingerprint matrix from SMILES list.

    优先使用 Aglycon_SMILES; 如果无效则回退到 canonical_smiles。
    设计原则: 我们聚类的是苷元, 不是整个分子, 因为同一骨架配不同糖链
    是我们要观察的现象 — 如果用全分子指纹, 糖链差异会干扰骨架聚类。

    Returns:
        (n_valid, nBits) 的 numpy 数组, 以及有效行索引列表
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem

    fps = []
    validIndices = []

    for i, smi in enumerate(smilesList):
        smi = str(smi).strip()
        # 移除 dummy atom 标记 (Remove dummy atom markers from cleavage)
        smi = smi.replace("[100*]", "").replace("[200*]", "")
        if not smi or smi == "nan":
            continue
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        try:
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nBits)
            arr = np.zeros(nBits, dtype=np.uint8)
            fp.SetBitsFromList(list(fp.GetOnBits()))
            # 转换为 numpy 数组 (Convert to array)
            for bit in fp.GetOnBits():
                arr[bit] = 1
            fps.append(arr)
            validIndices.append(i)
        except Exception:
            continue

    if not fps:
        return np.array([]), []

    return np.vstack(fps), validIndices


# =====================================================================
# 2. 数据准备与采样 (Data Preparation & Sampling)
# =====================================================================

def prepareData(df: pd.DataFrame, sampleSize: int = 10000) -> pd.DataFrame:
    """
    过滤、采样数据用于 UMAP 可视化。
    Filter and sample data for UMAP visualization.

    采样策略 (Sampling Strategy):
      - 保留所有稀有类 (如 Nucleotide Sugar, Glycopeptide)
      - 对大类做比例采样, 确保 Unclassified 不淹没小类
    """
    # 清理 Superclass 名称 (Clean Superclass)
    df["_class"] = df["Superclass"].apply(lambda x: cleanSuperclass(str(x)))

    # 过滤假苷元 (Filter out simple glycosides first)
    if "Is_Simple_Glycoside" in df.columns:
        nBefore = len(df)
        df = df[df["Is_Simple_Glycoside"].astype(str) != "True"].copy()
        print(f"  Filtered Is_Simple_Glycoside: {nBefore:,} → {len(df):,}")

    # 过滤无糖序列的异常值
    validDf = df[
        df["Sugar_Sequence"].notna()
        & (df["Sugar_Sequence"] != "")
        & (df["Sugar_Sequence"] != "nan")
    ].copy()

    # 过滤无苷元的行
    validDf = validDf[
        validDf["Aglycon_SMILES"].notna()
        & (validDf["Aglycon_SMILES"] != "")
        & (validDf["Aglycon_SMILES"] != "nan")
    ]

    print(f"  Valid rows (has Sugar + Aglycon): {len(validDf):,}")

    if len(validDf) <= sampleSize:
        return validDf.reset_index(drop=True)

    # 分层采样 (Stratified sampling)
    classCounts = validDf["_class"].value_counts()

    # 稀有类全量保留 (Keep all rare classes in full)
    rareThreshold = sampleSize // 20  # 500 for 10K sample
    rareClasses = classCounts[classCounts <= rareThreshold].index.tolist()
    rareDf = validDf[validDf["_class"].isin(rareClasses)]

    # 大类按比例采样 (Proportional sampling for large classes)
    commonDf = validDf[~validDf["_class"].isin(rareClasses)]
    remaining = sampleSize - len(rareDf)

    if remaining > 0 and len(commonDf) > 0:
        frac = remaining / len(commonDf)
        frac = min(1.0, frac)
        sampledCommon = commonDf.sample(n=min(remaining, len(commonDf)),
                                        random_state=42)
    else:
        sampledCommon = pd.DataFrame()

    result = pd.concat([rareDf, sampledCommon], ignore_index=True)
    print(f"  Sampled: {len(result):,} (rare={len(rareDf)}, common={len(sampledCommon)})")

    return result


def cleanSuperclass(val: str) -> str:
    """清理 Superclass (去 Tanimoto 后缀, 合并子类型)."""
    if not val or val == "nan":
        return "Unclassified"
    s = str(val).strip()
    if "(Tanimoto=" in s:
        s = s[:s.index("(Tanimoto=")].strip()
    # 合并 Glycolipid 子类型 (Merge Glycolipid subtypes)
    if s.startswith("Glycolipid"):
        s = "Glycolipids"
    return s if s else "Unclassified"


# =====================================================================
# 3. UMAP 降维 + 可视化 (UMAP Reduction + Visualization)
# =====================================================================

def runUmapAndPlot(
    sampleDf: pd.DataFrame,
    outputPath: str,
    nNeighbors: int = 30,
    minDist: float = 0.3,
    metric: str = "jaccard",
):
    """
    执行 UMAP 降维并生成交互式 Plotly 散点图。
    Run UMAP and generate interactive Plotly scatter plot.

    参数选择 (Parameter Rationale):
      - n_neighbors=30: 平衡局部/全局结构
      - min_dist=0.3: 允许一定散开, 避免过度挤压
      - metric=jaccard: 最适合二值指纹 (binary fingerprints)
    """
    import umap
    import plotly.express as px

    print("\n  [Step 1] Computing Morgan fingerprints...")
    t0 = time.time()

    # 优先用 Aglycon, 回退到 canonical_smiles
    smilesCol = "Aglycon_SMILES"
    fpMatrix, validIdx = computeMorganFingerprints(
        sampleDf[smilesCol].tolist()
    )

    if len(fpMatrix) == 0:
        print("  [ERROR] No valid fingerprints computed!")
        return

    # 重新对齐 DataFrame (Realign DataFrame to valid indices)
    plotDf = sampleDf.iloc[validIdx].copy().reset_index(drop=True)
    print(f"    FP matrix: {fpMatrix.shape} ({time.time()-t0:.1f}s)")

    # UMAP 降维 (UMAP dimensionality reduction)
    print(f"\n  [Step 2] UMAP embedding ({fpMatrix.shape[0]} points)...")
    t1 = time.time()
    reducer = umap.UMAP(
        n_neighbors=nNeighbors,
        min_dist=minDist,
        n_components=2,
        metric=metric,
        random_state=42,
        n_jobs=1,   # 避免 Windows multiprocessing 问题
        verbose=False,
    )
    embedding = reducer.fit_transform(fpMatrix)
    plotDf["UMAP_X"] = embedding[:, 0]
    plotDf["UMAP_Y"] = embedding[:, 1]
    print(f"    Done ({time.time()-t1:.1f}s)")

    # 准备糖链标签 (Prepare sugar labels: Top 10 + Other)
    topSugars = plotDf["Sugar_Sequence"].value_counts().head(10).index.tolist()
    plotDf["Sugar_Top10"] = plotDf["Sugar_Sequence"].apply(
        lambda s: s if s in topSugars else "Other"
    )

    # 清理 hover 数据 (Clean hover data)
    plotDf["Name_Short"] = plotDf["name"].fillna("").apply(
        lambda x: str(x)[:40] + "..." if len(str(x)) > 40 else str(x)
    )
    plotDf["Org_Short"] = plotDf["organisms"].fillna("").apply(
        lambda x: str(x).split("|")[0][:35] if x else ""
    )
    plotDf["Mods"] = plotDf["Glycan_Modifications"].fillna("None")

    # 配色方案 (Color scheme: Unclassified → grey, named classes → vivid)
    classOrder = plotDf["_class"].value_counts().index.tolist()
    # 将 Unclassified 放到最后 (Move Unclassified to end for rendering order)
    if "Unclassified" in classOrder:
        classOrder.remove("Unclassified")
        classOrder.append("Unclassified")

    colorMap = {}
    vividColors = [
        "#e74c3c", "#3498db", "#2ecc71", "#f39c12", "#9b59b6",
        "#1abc9c", "#e67e22", "#e91e63", "#00bcd4", "#8bc34a",
        "#ff5722", "#607d8b", "#795548", "#ffc107", "#673ab7",
        "#4caf50", "#ff9800", "#03a9f4", "#cddc39", "#f44336",
    ]
    for i, cls in enumerate(classOrder):
        if cls == "Unclassified":
            colorMap[cls] = "rgba(180, 180, 180, 0.3)"
        else:
            colorMap[cls] = vividColors[i % len(vividColors)]

    # Plotly 散点图 (Plotly scatter plot)
    print(f"\n  [Step 3] Generating Plotly scatter...")

    fig = px.scatter(
        plotDf,
        x="UMAP_X",
        y="UMAP_Y",
        color="_class",
        category_orders={"_class": classOrder},
        color_discrete_map=colorMap,
        hover_data={
            "Name_Short": True,
            "Sugar_Sequence": True,
            "Org_Short": True,
            "Mods": True,
            "_class": False,
            "UMAP_X": False,
            "UMAP_Y": False,
        },
        labels={
            "_class": "Superclass",
            "UMAP_X": "UMAP-1",
            "UMAP_Y": "UMAP-2",
            "Name_Short": "Name",
            "Org_Short": "Organism",
            "Mods": "Modifications",
        },
        title="Aglycon Chemical Space — UMAP of Morgan Fingerprints<br>"
              "<sub>Each dot = one glycoconjugate | Color = Superclass | "
              "Grey = Unclassified 'dark matter'</sub>",
        opacity=0.7,
    )

    fig.update_traces(
        marker=dict(size=4, line=dict(width=0)),
        selector=dict(mode="markers"),
    )
    # Unclassified 点更小更透明 (Make Unclassified smaller)
    fig.for_each_trace(
        lambda t: t.update(marker=dict(size=2, opacity=0.2))
        if t.name == "Unclassified" else ()
    )

    fig.update_layout(
        width=1400, height=900,
        font=dict(family="Segoe UI, Arial", size=12),
        paper_bgcolor="#fafafa",
        plot_bgcolor="#111122",
        legend=dict(
            title="Superclass",
            font=dict(size=10),
            itemsizing="constant",
        ),
        xaxis=dict(showgrid=False, zeroline=False),
        yaxis=dict(showgrid=False, zeroline=False),
    )

    fig.write_html(outputPath)
    print(f"  Saved: {outputPath}")

    # 统计信息 (Summary stats)
    print(f"\n  [Summary]")
    print(f"    Total points plotted: {len(plotDf):,}")
    for cls in classOrder[:10]:
        n = len(plotDf[plotDf["_class"] == cls])
        print(f"    {cls:30s} {n:,}")

    return plotDf  # 返回含 UMAP 坐标的 DataFrame


# =====================================================================
# Main
# =====================================================================

def main():
    parser = argparse.ArgumentParser(description="GlycoNP UMAP Chemical Space")
    parser.add_argument("--input", type=str, default=None, help="Input CSV")
    parser.add_argument("--sample", type=int, default=10000, help="Sample size")
    parser.add_argument("--neighbors", type=int, default=30, help="UMAP n_neighbors")
    parser.add_argument("--min-dist", type=float, default=0.3, help="UMAP min_dist")
    args = parser.parse_args()

    baseDir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    reportDir = os.path.join(baseDir, "reports")

    if args.input:
        inputPath = args.input
    else:
        # 优先使用清洗后数据 (Prefer cleaned data)
        inputPath = os.path.join(reportDir, "GlycoNP_Pipeline_Full_Cleaned.csv")
        if not os.path.exists(inputPath):
            inputPath = os.path.join(reportDir, "GlycoNP_Pipeline_Full.csv")

    print("=" * 70)
    print("  GlycoNP UMAP Chemical Space Clustering")
    print("=" * 70)
    print(f"  Input: {inputPath}")

    df = pd.read_csv(inputPath, low_memory=False, dtype=str, encoding="utf-8-sig")
    print(f"  Loaded {len(df):,} rows")

    sampleDf = prepareData(df, sampleSize=args.sample)

    # UMAP 1: Superclass 配色 (Original)
    outputPath = os.path.join(reportDir, "UMAP_Chemical_Space.html")
    umapDf = runUmapAndPlot(
        sampleDf, outputPath,
        nNeighbors=args.neighbors,
        minDist=args.min_dist,
    )

    # UMAP 2: ΔLogP 颜色映射 (ΔLogP color-mapped view)
    if umapDf is not None and "Delta_LogP" in umapDf.columns:
        plotDeltaLogPUmap(umapDf, reportDir)

    print(f"\n{'='*70}")
    print(f"  Done!")
    print(f"{'='*70}")


def plotDeltaLogPUmap(sampleDf: pd.DataFrame, reportDir: str):
    """
    生成按 ΔLogP 颜色映射的 UMAP 散点图。
    Generate ΔLogP color-mapped UMAP scatter plot.

    化学意义: 颜色越深(蓝) = 糖链极性改造越极端。
    """
    import plotly.express as px

    print(f"\n  [ΔLogP UMAP] Generating ΔLogP color-mapped view...")

    if "UMAP_X" not in sampleDf.columns:
        print("    [SKIP] No UMAP coordinates — run Superclass UMAP first")
        return

    plotDf = sampleDf.copy()
    plotDf["Delta_LogP_f"] = pd.to_numeric(
        plotDf["Delta_LogP"], errors="coerce")
    plotDf = plotDf[plotDf["Delta_LogP_f"].notna()].copy()

    if len(plotDf) == 0:
        print("    [SKIP] No valid ΔLogP values")
        return

    # 准备 hover (Prepare hover)
    plotDf["Name_Short"] = plotDf["name"].fillna("").apply(
        lambda x: str(x)[:40] + "..." if len(str(x)) > 40 else str(x)
    )

    fig = px.scatter(
        plotDf,
        x="UMAP_X", y="UMAP_Y",
        color="Delta_LogP_f",
        color_continuous_scale="RdYlBu_r",  # 红(低ΔLogP) → 蓝(高ΔLogP)
        range_color=[-5, 15],
        hover_data={
            "Name_Short": True,
            "Sugar_Sequence": True,
            "Delta_LogP_f": ":.2f",
            "UMAP_X": False,
            "UMAP_Y": False,
        },
        labels={
            "Delta_LogP_f": "ΔLogP",
            "UMAP_X": "UMAP-1",
            "UMAP_Y": "UMAP-2",
            "Name_Short": "Name",
        },
        title="Polarity Remodeling Landscape — ΔLogP on Aglycon UMAP<br>"
              "<sub>Red = low ΔLogP (sugar adds little polarity) | "
              "Blue = high ΔLogP (sugar dramatically increases water solubility)</sub>",
        opacity=0.7,
    )

    fig.update_traces(marker=dict(size=4, line=dict(width=0)))
    fig.update_layout(
        width=1400, height=900,
        font=dict(family="Segoe UI, Arial", size=12),
        paper_bgcolor="#fafafa",
        plot_bgcolor="#111122",
        xaxis=dict(showgrid=False, zeroline=False),
        yaxis=dict(showgrid=False, zeroline=False),
    )

    outPath = os.path.join(reportDir, "UMAP_DeltaLogP.html")
    fig.write_html(outPath)
    print(f"    Saved: {outPath}")
    print(f"    Points: {len(plotDf):,}")
    print(f"    Mean ΔLogP: {plotDf['Delta_LogP_f'].mean():.2f}")


if __name__ == "__main__":
    main()
