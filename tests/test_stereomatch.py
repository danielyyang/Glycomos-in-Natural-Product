"""
Test stereochemical matching using RDKit's HasSubstructMatch with useChirality=True.
"""
import sys, os
from rdkit import Chem
from rdkit.Chem import AllChem

SUGAR_SMILES_LIB = {
    "Glc": "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O", 
    "Gal": "OC[C@H]1O[C@@H](O)[C@@H](O)[C@@H](O)[C@@H]1O",
    "Man": "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O",
    "All": "OC[C@H]1O[C@@H](O)[C@H](O)[C@H](O)[C@@H]1O",   
}

maltose_smi = "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O[C@H]1OC(CO)[C@@H](O)[C@H](O)[C@H]1O"
mol = Chem.MolFromSmiles(maltose_smi)
Chem.AssignStereochemistry(mol, force=True)

# Important: To match stereochemistry, we should use identical rings.
# But SMARTS matching is better for sub-components (monosaccharides inside polymers).
# Let's convert our library to SMARTS query, or simply Mol objects acting as queries.

import copy

print(f"Testing Maltose: {maltose_smi}")
for name, smi in SUGAR_SMILES_LIB.items():
    query = Chem.MolFromSmiles(smi)
    # The problem is that the query has explicit OH groups, while the polymer has O-C
    # We must convert the query's OH groups to O-* (wildcards) or just O with any connections.
    
    # Let's adjust the query: any O atom should be allowed to have an extra connection
    q_sma = Chem.MolToSmarts(query)
    # It might be easier to build a smarts matcher manually, or replace OH with [O;X2] etc.
    # Actually, RDKit's MolFromSmarts doesn't keep stereochemistry easily without explicit notation.
    # What if we just replace all terminal H on O in the query molecule with Dummies?
    
    mod_query = Chem.MolFromSmiles(smi)
    for qatom in mod_query.GetAtoms():
        if qatom.GetSymbol() == 'O':
            # Allow oxygen to match any O in the target (whether OH or O-C)
            qatom.SetQuery(Chem.rdqueries.AtomNumEqualsQueryAtom(8))
            
    matches = mol.GetSubstructMatches(mod_query, useChirality=True)
    print(f"  Match {name}: {len(matches)} times")

# Wait, the R/S signature approach was actually very good and fast IF we could just look at the up/down bonds relative to the ring instead of formal CIP priorities!
# We can determine relative stereochemistry (up/down) ourselves based on geometry!
