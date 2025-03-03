import pandas as pd

# Path to cleaned dataset
file_path = r"C:\Users\mh3511\Desktop\drug_discovery_ai\data\cleaned_data.csv"

# Load and display columns
df = pd.read_csv(file_path)
print("Columns in dataset:", df.columns)
print(df.head())
