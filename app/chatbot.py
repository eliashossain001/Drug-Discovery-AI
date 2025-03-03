import openai
import os
from fastapi import FastAPI
from pydantic import BaseModel

# Initialize FastAPI chatbot app
app = FastAPI(title="Drug Discovery Chatbot")

# Load OpenAI API Key from environment variable
openai.api_key = os.getenv("OPENAI_API_KEY")

# Define request model
class ChatRequest(BaseModel):
    question: str

# Function to query OpenAI's GPT-4 or GPT-3.5 Turbo
def get_gpt_response(question):
    """
    Queries OpenAI's GPT-4 (or GPT-3.5 Turbo) for answering drug-related questions.
    """
    try:
        response = openai.ChatCompletion.create(
            model="gpt-4",  # Change to "gpt-3.5-turbo" if preferred
            messages=[
                {"role": "system", "content": "You are an AI assistant specializing in drug discovery and QSAR analysis."},
                {"role": "user", "content": question}
            ]
        )
        return response["choices"][0]["message"]["content"]
    except Exception as e:
        return f"Error querying OpenAI: {e}"

# Chatbot API endpoint
@app.post("/chat/")
def chat_with_gpt(request: ChatRequest):
    """
    Handles chat queries and returns AI-generated responses.
    """
    answer = get_gpt_response(request.question)
    return {"question": request.question, "answer": answer}

# Root endpoint
@app.get("/")
def root():
    return {"message": "Welcome to the Drug Discovery Chatbot API. Use /chat to interact."}
