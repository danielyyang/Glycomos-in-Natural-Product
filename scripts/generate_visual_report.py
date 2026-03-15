"""
Phase 7 Test: Generate visual report for 100 sample compounds
Full pipeline: Phase 2 (cleavage) -> Phase 5 (features) -> Phase 6 (classify) -> Phase 7 (visualize)
"""
import os, sys, time
import pandas as pd
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from rdkit import Chem
from lib.sugar_utils import find_mapped_sugar_units
from lib.glycosidic_cleavage import cleaveWithConservation
from lib.phase5_features import characterizeGlycan, extractMurckoScaffold
from lib.phase6_classifier import classifyAglycon, buildReferenceFingerprints
from lib.phase7_visualizer import generateHtmlReport


def main():
    baseDir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    dataPath = os.path.join(baseDir, "reports", "Coconut_Sugar_Check.csv")
    outputHtml = os.path.join(baseDir, "reports", "Sample_Visual_Report.html")

    print("=" * 70)
    print("  Phase 7: Visual Report Generation — 100 Compound Pipeline")
    print("=" * 70)

    df = pd.read_csv(dataPath, low_memory=False, dtype=str, encoding="utf-8-sig")
    print(f"  Total rows: {len(df)}")

    # Sample 100 valid compounds
    validMask = df["canonical_smiles"].notna() & (~df["canonical_smiles"].isin(["", "nan", "NULL"]))
    sampleDf = df[validMask].sample(n=100, random_state=7).copy()
    print(f"  Sampled: {len(sampleDf)} compounds")

    # Build Tanimoto reference from the broader dataset
    print("\n[Step 1] Building Tanimoto reference...")
    t0 = time.time()
    classCol = "np_classifier_superclass"
    existing = df[validMask & df[classCol].notna() & (~df[classCol].isin(["", "nan", "NULL"]))].sample(
        n=min(500, len(df)), random_state=42
    )
    refFps, refLabels = buildReferenceFingerprints(existing, "canonical_smiles", classCol)
    print(f"  Reference: {len(refFps)} FPs in {time.time()-t0:.1f}s")

    # Run full pipeline on each compound
    print("\n[Step 2] Running Phase 2 -> 5 -> 6 pipeline...")
    t1 = time.time()

    newCols = {
        "sugar_sequence": [], "sugar_functional_group": [],
        "murcko_scaffold": [], "classification": [], "glycolipid_flag": [],
    }

    for i, (_, row) in enumerate(sampleDf.iterrows()):
        smiles = str(row["canonical_smiles"])
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                raise ValueError("Invalid SMILES")

            # Phase 2: cleavage
            units = find_mapped_sugar_units(mol)
            glycan, aglycon, meta = cleaveWithConservation(mol, units)

            # Phase 5: glycan sequence + aglycon scaffold
            glycanFeatures = characterizeGlycan(glycan)
            scaffold = extractMurckoScaffold(aglycon)

            # Phase 6: classification
            classResult = classifyAglycon(smiles, refFps, refLabels)

            newCols["sugar_sequence"].append(glycanFeatures.get("sugar_sequence", ""))
            newCols["sugar_functional_group"].append(glycanFeatures.get("sugar_functional_group", ""))
            newCols["murcko_scaffold"].append(scaffold)
            newCols["classification"].append(classResult.get("classification", ""))
            newCols["glycolipid_flag"].append(classResult.get("glycolipid_flag", ""))

        except Exception as e:
            newCols["sugar_sequence"].append("")
            newCols["sugar_functional_group"].append("")
            newCols["murcko_scaffold"].append("")
            newCols["classification"].append(f"Error: {e}")
            newCols["glycolipid_flag"].append("")

        if (i + 1) % 25 == 0:
            print(f"    Processed {i+1}/100")

    for col, values in newCols.items():
        sampleDf[col] = values

    pipelineTime = time.time() - t1
    print(f"  Pipeline done in {pipelineTime:.1f}s ({pipelineTime/100*1000:.0f}ms/compound)")

    # Phase 7: generate HTML report
    print("\n[Step 3] Generating HTML report with highlighted structures...")
    t2 = time.time()
    generateHtmlReport(
        sampleDf, outputHtml,
        smilesCol="canonical_smiles",
        inchikeyCol="standard_inchi_key",
        sugarSeqCol="sugar_sequence",
        scaffoldCol="murcko_scaffold",
        classCol="classification",
        glycolipidCol="glycolipid_flag",
        maxRows=100,
        imgSize=(350, 250),
    )
    reportTime = time.time() - t2

    # Summary stats
    print(f"\n{'='*70}")
    print("  RESULTS")
    print(f"{'='*70}")

    classDist = sampleDf["classification"].value_counts().head(8)
    print(f"\n  [Class Distribution (top 8)]")
    for cls, count in classDist.items():
        print(f"    {str(cls):45} {count:3}")

    glycoCount = sampleDf[sampleDf["glycolipid_flag"] != ""].shape[0]
    seqCount = sampleDf[sampleDf["sugar_sequence"] != ""].shape[0]
    print(f"\n  Sugar sequences generated: {seqCount}/100")
    print(f"  Glycolipids flagged: {glycoCount}/100")
    print(f"  Report time: {reportTime:.1f}s")
    print(f"  Total time: {time.time()-t0:.1f}s")
    print(f"  Output: {outputHtml}")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
