from rdkit import Chem
import pandas as pd
import glob

# Find the Excel file that starts with "PO-"
excel_files = glob.glob('PO-*.xlsx')
if not excel_files:
    raise FileNotFoundError("No Excel file starting with 'PO-' found.")
excel_file = excel_files[0]  # Use the first matching file

# 'Comment'&'Amount' are the colums from PO-Excel file I want to merge to my Output.sdf
df_excel = pd.read_excel(excel_file, usecols=['ID', 'Comment','Amount'])

# Convert ID in Excel file to string
df_excel['ID'] = df_excel['ID'].astype(str)

# Load the SDF file of your quote to Enamine (usually contains more cmpds than in PO-Excel file)
sdf_file = 'Your_quote.sdf'
sdf_supplier = Chem.SDMolSupplier(sdf_file, sanitize=False, removeHs=False)

# Output file
output_file = 'Merged_for_biology.sdf'
writer = Chem.SDWriter(output_file)

corrected_compounds = []

for mol in sdf_supplier:
    if mol is not None:
        bri_id = mol.GetProp('BRI_ID')
        bri_id = str(bri_id)

        try:
            # Attempt to sanitize the molecule
            Chem.SanitizeMol(mol)
        except:
            # Record the ID of the molecule that needed correction
            corrected_compounds.append(bri_id)

        # Check if BRI_ID is in the PO-Excel file
        if bri_id in df_excel['ID'].values:
            # Get molecule property of your interest to added to your output.sdf
            comment = str(df_excel.loc[df_excel['ID'] == bri_id, 'Comment'].values[0])
            amount = str(df_excel.loc[df_excel['ID'] == bri_id, 'Amount'].values[0])
            mol.SetProp('Comment', comment)
            mol.SetProp('Amount', amount)

        writer.write(mol)

writer.close()

print("SDF file has been updated and saved.")
print("Compounds with corrected valences:", corrected_compounds)