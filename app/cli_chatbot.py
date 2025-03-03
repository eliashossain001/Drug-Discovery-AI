import requests

API_URL = "http://127.0.0.1:8000/chat/"

def ask_question():
    """
    Allows users to ask drug-related questions via the command line.
    """
    print("💬 Drug Discovery Chatbot (Type 'exit' to quit)\n")
    
    while True:
        user_input = input("You: ")
        if user_input.lower() == "exit":
            print("Goodbye! 👋")
            break

        response = requests.post(API_URL, json={"question": user_input})
        if response.status_code == 200:
            answer = response.json().get("answer", "Sorry, I couldn't process that.")
            print(f"Bot: {answer}\n")
        else:
            print("Error: Could not get a response from the chatbot.")

if __name__ == "__main__":
    ask_question()
