import sys, os
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
sys.path.append(os.path.abspath(os.path.dirname(os.path.dirname(__file__))))
from lib import sugar_utils, sugar_sequence
from tests.test_linkage import TEST_CASES

os.makedirs('tests/images', exist_ok=True)

for name, smi in TEST_CASES.items():
    mol = Chem.MolFromSmiles(smi)
    if not mol: 
        print(f"Skipping {name}: Invalid SMILES")
        continue
    
    # 1. Generate units
    try:
        units, atom_map = sugar_utils.get_sugar_units(mol)
        seq, mods = sugar_sequence.generate_refined_sequence(mol)
    except Exception as e:
        seq = f"Error: {e}"
        units = []
        
    highlight_atoms = []
    highlight_colors = {}
    
    # 2. Add properties & highlights
    for u in units:
        # Paint rings Light Blue
        for idx in u['ring_atoms']:
            highlight_atoms.append(idx)
            highlight_colors[idx] = (0.5, 0.8, 1.0)
            
        o_idx = u['ring_oxygen']
        c1_idx = u['anomeric_idx']
        
        # We will label the C1 atom clearly with the recognized sugar name!
        c1_atom = mol.GetAtomWithIdx(c1_idx)
        label = f"{u['name']}_{u['id']} ({u['anomeric_config']})"
        c1_atom.SetProp('atomNote', label)
        
        # Color C1 orange to make it easy to see
        highlight_colors[c1_idx] = (1.0, 0.6, 0.0)
        
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
        
    # 3. Draw image
    d2d = rdMolDraw2D.MolDraw2DCairo(1200, 800)
    dopts = d2d.drawOptions()
    dopts.addAtomIndices = True # Display atom numbers for debugging connectivity
    dopts.useBWAtomPalette()
    dopts.padding = 0.1
    dopts.legendFontSize = 24
    
    d2d.DrawMolecule(mol, legend=f"{name}\nSequence: {seq}", highlightAtoms=highlight_atoms, highlightAtomColors=highlight_colors)
    d2d.FinishDrawing()
    
    safe_name = name.split(' ')[0].replace('(', '').replace(')', '')
    img_path = f"tests/images/{safe_name}.png"
    with open(img_path, "wb") as f:
        f.write(d2d.GetDrawingText())
    print(f"Generated {img_path}")
