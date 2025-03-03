import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors

class SMILESFeatureExtractor:
    def __init__(self, input_file, output_file):
        """
        Initializes the SMILES feature extractor with input and output file paths.
        """
        self.input_file = input_file
        self.output_file = output_file

    def load_data(self):
        """
        Loads the cleaned dataset containing SMILES notation and pIC50 values.
        """
        try:
            df = pd.read_csv(self.input_file)
            if "smiles" not in df.columns:
                raise ValueError("SMILES column not found in the dataset.")
            if "pic50" not in df.columns:
                raise ValueError("Target column 'pIC50' not found in the dataset.")
            print(f"Data loaded successfully! Shape: {df.shape}")
            return df
        except Exception as e:
            print(f"Error loading data: {e}")
            return None

    def extract_features(self, df):
        """
        Extracts molecular features from the SMILES notation and keeps pIC50 values.
        """
        features = []
        for index, row in df.iterrows():
            smiles = row["smiles"]
            pIC50 = row["pic50"]  # Keep pIC50 value
            mol = Chem.MolFromSmiles(smiles)

            if mol:
                mol_weight = Descriptors.MolWt(mol)
                num_h_donors = Descriptors.NumHDonors(mol)
                num_h_acceptors = Descriptors.NumHAcceptors(mol)
                num_rotatable_bonds = Descriptors.NumRotatableBonds(mol)
                logp = Descriptors.MolLogP(mol)  # Lipophilicity (LogP)
                
                features.append([smiles, mol_weight, num_h_donors, num_h_acceptors, num_rotatable_bonds, logp, pIC50])
            else:
                print(f"Invalid SMILES string at index {index}: {smiles}")

        # Create DataFrame with extracted features
        feature_df = pd.DataFrame(features, columns=[
            "smiles", "molecular_weight", "num_h_donors", "num_h_acceptors", "num_rotatable_bonds", "logp", "pIC50"
        ])

        print(f"Extracted features shape: {feature_df.shape}")
        return feature_df

    def save_features(self, feature_df):
        """
        Saves the extracted molecular features to a CSV file.
        """
        try:
            feature_df.to_csv(self.output_file, index=False)
            print(f"Feature-engineered data saved to {self.output_file}")
        except Exception as e:
            print(f"Error saving feature data: {e}")

    def run(self):
        """
        Executes the full feature extraction pipeline.
        """
        df = self.load_data()
        if df is not None:
            feature_df = self.extract_features(df)
            self.save_features(feature_df)

# Execute script
if __name__ == "__main__":
    input_path = r"C:\Users\mh3511\Desktop\drug_discovery_ai\data\cleaned_data.csv"  # Path to cleaned dataset
    output_path = r"C:\Users\mh3511\Desktop\drug_discovery_ai\data\feature_engineered_data.csv"

    extractor = SMILESFeatureExtractor(input_path, output_path)
    extractor.run()
