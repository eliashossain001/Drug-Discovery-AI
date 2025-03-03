import pandas as pd
import json

# Load dataset
df = pd.read_csv(r"C:\Users\mh3511\Desktop\drug_discovery_ai\data\cleaned_data.csv")

# Print available columns
print("Available columns:", df.columns)

# Adjust these column names based on what the debug output shows
QUESTION_COLUMN = "your_actual_question_column_name"  # Change this!
ANSWER_COLUMN = "your_actual_answer_column_name"  # Change this!

# Prepare messages for fine-tuning
training_data = []

for _, row in df.iterrows():
    user_question = row[QUESTION_COLUMN]
    bot_response = row[ANSWER_COLUMN]

    training_data.append({
        "messages": [
            {"role": "system", "content": "You are MedXpert, an AI expert in drug discovery."},
            {"role": "user", "content": user_question},
            {"role": "assistant", "content": bot_response}
        ]
    })

# Save as JSONL format
output_file = "data/fine_tune_data.jsonl"
with open(output_file, "w") as f:
    for entry in training_data:
        f.write(json.dumps(entry) + "\n")

print(f"✅ Fine-tuning dataset created: {output_file}")
