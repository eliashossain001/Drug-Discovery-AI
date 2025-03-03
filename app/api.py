from fastapi import FastAPI
from pydantic import BaseModel
import joblib
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors

# Initialize FastAPI app
app = FastAPI(title="QSAR Drug Discovery API", version="1.0", description="API to predict pIC50 for drug candidates")

# Load the trained QSAR model
model_path = r"C:\Users\mh3511\Desktop\drug_discovery_ai\models\qsar_model.pkl"  # Adjust path if needed
model = joblib.load(model_path)

# Define request models
class SingleDrugRequest(BaseModel):
    smiles: str

class CompareDrugsRequest(BaseModel):
    smiles1: str
    smiles2: str

# Feature extraction function
def extract_features(smiles):
    """
    Converts a SMILES string into molecular descriptors.
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
    return None

# Prediction function
def predict_pic50(smiles):
    """
    Predicts the pIC50 value for a given SMILES string.
    """
    features = extract_features(smiles)
    if features:
        features_df = pd.DataFrame([features], columns=[
            "molecular_weight", "num_h_donors", "num_h_acceptors", "num_rotatable_bonds", "logp"
        ])
        prediction = model.predict(features_df)[0]
        return round(prediction, 4)
    return None

# Single drug prediction endpoint
@app.post("/predict/", tags=["Prediction"])
def predict_drug(data: SingleDrugRequest):
    """
    Predicts the pIC50 for a given drug (SMILES input).
    """
    pic50 = predict_pic50(data.smiles)
    if pic50 is not None:
        return {"smiles": data.smiles, "predicted_pIC50": pic50}
    return {"error": "Invalid SMILES string"}

# Compare two drugs endpoint
@app.post("/compare/", tags=["Comparison"])
def compare_drugs(data: CompareDrugsRequest):
    """
    Compares two drug candidates and returns the stronger one based on predicted pIC50.
    """
    pic50_1 = predict_pic50(data.smiles1)
    pic50_2 = predict_pic50(data.smiles2)

    if pic50_1 is None or pic50_2 is None:
        return {"error": "Invalid SMILES string provided"}

    result = {
        "drug_1": {"smiles": data.smiles1, "predicted_pIC50": pic50_1},
        "drug_2": {"smiles": data.smiles2, "predicted_pIC50": pic50_2},
        "stronger_candidate": "Drug 1" if pic50_1 > pic50_2 else "Drug 2" if pic50_2 > pic50_1 else "Both drugs are equally effective"
    }
    return result

# Root endpoint
@app.get("/", tags=["General"])
def root():
    """
    Root API endpoint.
    """
    return {"message": "Welcome to the QSAR Drug Discovery API. Use /docs for API testing."}
