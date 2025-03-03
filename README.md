# 💊 MedXpert - AI-Powered Drug Discovery & Comparison  

MedXpert is an **AI-driven chatbot** designed to assist researchers, pharmacists, and healthcare professionals in **drug discovery, QSAR modeling, and drug comparison**. It leverages **OpenAI’s GPT-4 API** to provide **real-time insights** into drug interactions, mechanisms, and pharmacology. Additionally, the **Drug Comparison module** allows users to analyze two drugs based on their molecular properties and predict efficacy.  

---

## 🚀 **Features**
✅ **MedXpert Chatbot** - AI-powered chatbot for answering pharma-related queries.  
✅ **Drug Comparison API** - Compares two drugs using QSAR-based analysis.  
✅ **Interactive Web UI** - Chatbot with an intuitive user interface.  
✅ **FastAPI Backend** - High-performance API for drug-related queries.  
✅ **SMILES-Based Analysis** - Converts molecular structures into numerical data.  
✅ **Locally Testable** - Supports **local API requests & command-line interaction**.  

---

## 📂 **Project Structure**
drug_discovery_ai/ │── app/ │ ├── api.py # FastAPI API logic │ ├── chat_ui.py # Web-based chatbot UI │ 
├── chatbot.py # Chatbot backend logic │ ├── cli_chatbot.py # CLI chatbot for terminal use │── data/ │ 
├── cleaned_data.csv # Preprocessed dataset │ ├── feature_engineered_data.csv # Features for drug comparison │── models/ │ 
├── qsar_model.pkl # Pretrained QSAR model │── src/ │ ├── inference/ │ │ ├── compare_drugs.py # Drug comparison logic │ │
├── predict_activity.py # QSAR model predictions │ ├── preprocessing/ │ │ ├── smiles_to_features.py # Converts SMILES notation to features │
├── training/ │ │ ├── train_qsar.py # Model training script │── utils/ │ ├── config.py # Centralized configuration file │
├── file_manager.py # Utility functions for file handling │── .env # API keys & environment variables
│── requirements.txt # Python dependencies │── README.md # Project documentation


---

## 🛠️ **Installation**
### **Clone the Repository**
```
git clone https://github.com/yourusername/drug_discovery_ai.git
cd drug_discovery_ai
```

## 🛠️ **Create Virtual Environment**
```
python -m venv venv
source venv/bin/activate  # For Mac/Linux
venv\Scripts\activate     # For Windows
```

## 🛠️ **Install Dependencies**
```
pip install -r requirements.txt
```

## 🛠️ **SET UP API KEY**
```
OPENAI_API_KEY="your-api-key-here"
```

## 🛠️ **How to Use?**
1️⃣ un the chatbot API:
```
uvicorn app.chat_ui:app --reload
```

2️⃣ Open in a browser:
👉 http://127.0.0.1:8000/chat-ui

3️⃣ Start chatting with MedXpert! 🤖💊

## 🛠️ **Drug Comparision API**
```
uvicorn app.api:app --reload
```
## Example API Request (Compare Two Drugs)
```
curl -X 'POST' \
  'http://127.0.0.1:8000/compare/' \
  -H 'Content-Type: application/json' \
  -d '{
    "smiles1": "CC(=O)Oc1ccccc1C(=O)O",
    "smiles2": "CC(C)Cc1ccc(cc1)C(C)C(=O)O"
}'
```

✅ The system will predict which drug has higher efficacy!

## Command-Line Interface (CLI) Chatbot
```
python app/cli_chatbot.py
```
Start asking questions like:

```
> What is QSAR modeling?
> How does LogP affect drug absorption?
```
## Configuration (config.py)
The config.py file manages project settings like:

* API keys
* Model paths
* File directories

## Example Usages
```
from utils.config import config

print("Using API Key:", config.OPENAI_API_KEY)
```

## 🎯 Contributing
Want to improve MedXpert? Fork the repo & create a pull request!
For major changes, open an issue first to discuss your ideas.

