import pandas as pd
from rdkit import Chem

file_path = 'PO-.xlsx'  # Replace with Emanine recieved file

df = pd.read_excel(file_path)

# Function to convert SMILES to canonical format
def convert_smiles(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        return Chem.MolToSmiles(mol, isomericSmiles=False, canonical=True, allHsExplicit=False)
    except:
        return None

# Convert and update SMILES strings in the DF
df['SMILES'] = df['SMILES'].apply(convert_smiles)

# Save the updated DF to a new Excel file
output_file_path = 'PO_SMILES_converted.xlsx'  
df.to_excel(output_file_path, index=False)

print(f"File saved to {output_file_path}")