import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import PandasTools, rdFMCS, AllChem

def create_fingerprint(mol):
    return AllChem.GetMorganFingerprint(mol, 2)

# Load the SDF files into DataFrames
df_converted_po = PandasTools.LoadSDF('converted_PO.sdf')
df_enamine_order = PandasTools.LoadSDF('Your_quote.sdf')

# Create fingerprints for each molecule
df_converted_po['Fingerprint'] = df_converted_po['ROMol'].apply(create_fingerprint)
df_enamine_order['Fingerprint'] = df_enamine_order['ROMol'].apply(create_fingerprint)

def structures_match(mol1, mol2, fp1, fp2):
    # First check using fingerprints
    similarity = DataStructs.TanimotoSimilarity(fp1, fp2)
    if similarity < 0.9: #sometimes NO2, or other functionalities are defined with the wrong valency in Enamine file
        return False
    # Perform MCS search if fingerprints are similar
    mcs = rdFMCS.FindMCS([mol1, mol2], timeout=10)  
    return mcs.numAtoms >= min(mol1.GetNumAtoms(), mol2.GetNumAtoms())

for i, row_po in df_converted_po.iterrows():
    for j, row_enamine in df_enamine_order.iterrows():
        if structures_match(row_po['ROMol'], row_enamine['ROMol'], row_po['Fingerprint'], row_enamine['Fingerprint']):
            # Combine data from both rows
            for col in df_enamine_order.columns:
                if col not in df_converted_po.columns:
                    df_converted_po.loc[i, col] = row_enamine[col]

# Save the merged DataFrame to a new SDF file
PandasTools.WriteSDF(df_converted_po, 'merged.sdf', molColName='ROMol', properties=list(df_converted_po.columns))
