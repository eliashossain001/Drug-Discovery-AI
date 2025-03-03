import pandas as pd
import os
import joblib
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error
from math import sqrt  # Import sqrt function for RMSE calculation

class QSARModelTrainer:
    def __init__(self, input_file, model_output):
        """
        Initializes the QSAR model trainer with input dataset and output model path.
        """
        self.input_file = input_file
        self.model_output = model_output

    def load_data(self):
        """
        Loads the feature-engineered dataset.
        """
        try:
            df = pd.read_csv(self.input_file)
            if "pIC50" not in df.columns:
                raise ValueError("Target column 'pIC50' not found in dataset.")
            print(f"Data loaded successfully! Shape: {df.shape}")
            return df
        except Exception as e:
            print(f"Error loading data: {e}")
            return None

    def prepare_data(self, df):
        """
        Splits data into features (X) and target variable (y), and performs train-test split.
        """
        # Define feature columns
        feature_columns = ["molecular_weight", "num_h_donors", "num_h_acceptors", "num_rotatable_bonds", "logp"]

        if not all(col in df.columns for col in feature_columns):
            raise ValueError("One or more required feature columns are missing from the dataset.")

        X = df[feature_columns]  # Features
        y = df["pIC50"]  # Target variable (bioactivity)

        # Split data (80% train, 20% test)
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

        print(f"Training set: {X_train.shape}, Test set: {X_test.shape}")
        return X_train, X_test, y_train, y_test

    def train_model(self, X_train, y_train):
        """
        Trains a RandomForestRegressor model.
        """
        model = RandomForestRegressor(n_estimators=100, random_state=42)
        model.fit(X_train, y_train)
        print("Model training complete.")
        return model

    def evaluate_model(self, model, X_test, y_test):
        """
        Evaluates the model using RMSE.
        """
        predictions = model.predict(X_test)
        mse = mean_squared_error(y_test, predictions)  # Compute Mean Squared Error (MSE)
        rmse = sqrt(mse)  # Compute Root Mean Squared Error (RMSE)
        print(f"Model RMSE: {rmse:.4f}")
        return rmse

    def save_model(self, model):
        """
        Saves the trained model to a file.
        """
        try:
            joblib.dump(model, self.model_output)
            print(f"Model saved to {self.model_output}")
        except Exception as e:
            print(f"Error saving model: {e}")

    def run(self):
        """
        Executes the full model training pipeline.
        """
        df = self.load_data()
        if df is not None:
            X_train, X_test, y_train, y_test = self.prepare_data(df)
            model = self.train_model(X_train, y_train)
            self.evaluate_model(model, X_test, y_test)
            self.save_model(model)

# Execute script
if __name__ == "__main__":
    input_path = r"C:\Users\mh3511\Desktop\drug_discovery_ai\data\feature_engineered_data.csv"  # Path to processed dataset
    model_path = r"C:\Users\mh3511\Desktop\drug_discovery_ai\models\qsar_model.pkl"  # Save trained model

    trainer = QSARModelTrainer(input_path, model_path)
    trainer.run()
