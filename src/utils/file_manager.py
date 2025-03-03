import os
import pandas as pd
import pickle

class FileManager:
    @staticmethod
    def load_csv(file_path):
        """Loads a CSV file into a Pandas DataFrame."""
        if os.path.exists(file_path):
            return pd.read_csv(file_path)
        else:
            raise FileNotFoundError(f"File not found: {file_path}")

    @staticmethod
    def save_pickle(obj, file_path):
        """Saves an object as a pickle file."""
        with open(file_path, "wb") as f:
            pickle.dump(obj, f)

    @staticmethod
    def load_pickle(file_path):
        """Loads a pickle file."""
        if os.path.exists(file_path):
            with open(file_path, "rb") as f:
                return pickle.load(f)
        else:
            raise FileNotFoundError(f"Pickle file not found: {file_path}")

file_manager = FileManager()
