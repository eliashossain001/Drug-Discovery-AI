import pandas as pd
import os

class DataCleaner:
    def __init__(self, input_file, output_file):
        """
        Initializes the DataCleaner with input and output file paths.
        """
        self.input_file = input_file
        self.output_file = output_file

    def load_data(self):
        """
        Loads the dataset into a Pandas DataFrame.
        """
        try:
            df = pd.read_csv(self.input_file)
            print(f"Data loaded successfully! Shape: {df.shape}")
            return df
        except Exception as e:
            print(f"Error loading data: {e}")
            return None

    def clean_data(self, df):
        """
        Cleans the dataset by handling missing values and standardizing column names.
        """
        # Standardize column names (convert to lowercase and replace spaces with underscores)
        df.columns = df.columns.str.lower().str.replace(" ", "_")

        # Ensure `pIC50` is numerical and remove invalid values
        if "pic50" in df.columns:
            df = df[df["pic50"] != "BLINDED"]  # Remove rows with "BLINDED"
            df["pic50"] = pd.to_numeric(df["pic50"], errors="coerce")  # Convert to numeric
            df = df.dropna(subset=["pic50"])  # Drop rows with NaN values in pIC50

        print(f"Cleaned data shape: {df.shape}")
        return df

    def save_cleaned_data(self, df):
        """
        Saves the cleaned dataset to the output file.
        """
        try:
            df.to_csv(self.output_file, index=False)
            print(f"Cleaned data saved to {self.output_file}")
        except Exception as e:
            print(f"Error saving cleaned data: {e}")

    def run(self):
        """
        Executes the full data cleaning pipeline.
        """
        df = self.load_data()
        if df is not None:
            df = self.clean_data(df)
            self.save_cleaned_data(df)

# Execute script
if __name__ == "__main__":
    input_path = r"C:\Users\mh3511\Desktop\drug_discovery_ai\data\DDH Data with Properties.csv"  # Adjust if needed
    output_path = r"C:\Users\mh3511\Desktop\drug_discovery_ai\data\cleaned_data.csv"

    cleaner = DataCleaner(input_path, output_path)
    cleaner.run()
