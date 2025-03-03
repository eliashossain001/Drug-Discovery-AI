import pandas as pd
import joblib
import sys
from rdkit import Chem
from rdkit.Chem import Descriptors

class DrugActivityPredictor:
    def __init__(self, model_path):
        """
        Initializes the predictor by loading the trained QSAR model.
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

# Execute script
if __name__ == "__main__":
    model_path = r"C:\Users\mh3511\Desktop\drug_discovery_ai\models\qsar_model.pkl"  # Path to the trained model
    predictor = DrugActivityPredictor(model_path)

    # Example usage (replace with actual SMILES string)
    test_smiles = "CC(=O)Oc1ccccc1C(=O)O"  # Example: Aspirin
    predicted_pIC50 = predictor.predict_activity(test_smiles)

    if predicted_pIC50 is not None:
        print(f"Predicted pIC50 for {test_smiles}: {predicted_pIC50}")
    else:
        print("Failed to predict activity.")
