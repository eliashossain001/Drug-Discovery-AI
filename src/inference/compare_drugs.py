import pandas as pd
import joblib
import sys
from rdkit import Chem
from rdkit.Chem import Descriptors

class DrugComparator:
    def __init__(self, model_path):
        """
        Initializes the comparator by loading the trained QSAR model.
        """
        self.model_path = model_path
        self.model = self.load_model()

    def load_model(self):
        """
        Loads the trained QSAR model.
        """
        try:
            model = joblib.load(self.model_path)
            print("Model loaded successfully.")
            return model
        except Exception as e:
            print(f"Error loading model: {e}")
            sys.exit(1)

    def extract_features_from_smiles(self, smiles):
        """
        Converts a SMILES string to molecular descriptors.
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            return [
                Descriptors.MolWt(mol),                # Molecular Weight
                Descriptors.NumHDonors(mol),           # Hydrogen Bond Donors
                Descriptors.NumHAcceptors(mol),        # Hydrogen Bond Acceptors
                Descriptors.NumRotatableBonds(mol),    # Rotatable Bonds
                Descriptors.MolLogP(mol)               # Lipophilicity (LogP)
            ]
        else:
            print(f"Invalid SMILES string: {smiles}")
            return None

    def predict_activity(self, smiles):
        """
        Predicts the pIC50 value for a given SMILES string.
        """
        features = self.extract_features_from_smiles(smiles)
        if features:
            features_df = pd.DataFrame([features], columns=[
                "molecular_weight", "num_h_donors", "num_h_acceptors", "num_rotatable_bonds", "logp"
            ])
            prediction = self.model.predict(features_df)[0]  # Predict pIC50
            return round(prediction, 4)  # Return rounded prediction
        return None

    def compare_drugs(self, smiles1, smiles2):
        """
        Compares two drug candidates based on their predicted pIC50 values.
        """
        pic50_1 = self.predict_activity(smiles1)
        pic50_2 = self.predict_activity(smiles2)

        if pic50_1 is None or pic50_2 is None:
            print("Error: One or both SMILES strings are invalid.")
            return

        print(f"Predicted pIC50 for Drug 1 ({smiles1}): {pic50_1}")
        print(f"Predicted pIC50 for Drug 2 ({smiles2}): {pic50_2}")

        if pic50_1 > pic50_2:
            print("✅ Drug 1 is the stronger antiviral candidate.")
        elif pic50_2 > pic50_1:
            print("✅ Drug 2 is the stronger antiviral candidate.")
        else:
            print("Both drugs have equal predicted bioactivity.")

# Execute script
if __name__ == "__main__":
    model_path = r"C:\Users\mh3511\Desktop\drug_discovery_ai\models\qsar_model.pkl"  # Path to trained QSAR model
    comparator = DrugComparator(model_path)

    # Example SMILES inputs (Aspirin vs. Ibuprofen)
    smiles_1 = "CC(=O)Oc1ccccc1C(=O)O"  # Aspirin
    smiles_2 = "CC(C)Cc1ccc(cc1)C(C)C(=O)O"  # Ibuprofen

    comparator.compare_drugs(smiles_1, smiles_2)
