from fastapi import FastAPI, Request
from fastapi.responses import HTMLResponse, JSONResponse
import openai
import os

# Initialize FastAPI app
app = FastAPI(title="Drug Discovery Chatbot")

# 🔒 Use API Key from Environment Variable
openai.api_key = os.getenv("OPENAI_API_KEY")

# Improved HTML page for better chat UI
html_content = """
<!DOCTYPE html>
<html>
<head>
    <title>💊 Drug Discovery Chatbot</title>
    <style>
        body { font-family: Arial, sans-serif; background-color: #f4f4f4; text-align: center; }
        .chat-container { width: 60%; margin: auto; background: white; padding: 20px; border-radius: 10px; box-shadow: 0px 0px 10px rgba(0, 0, 0, 0.1); }
        .chat-box { height: 500px; overflow-y: auto; border: 1px solid #ccc; padding: 15px; background: #fafafa; text-align: left; border-radius: 10px; display: flex; flex-direction: column; }
        .user-message, .bot-message { padding: 10px; border-radius: 15px; display: inline-block; max-width: 70%; }
        .user-message { background: #007bff; color: white; align-self: flex-end; }
        .bot-message { background: #e9ecef; color: black; align-self: flex-start; }
        .input-container { margin-top: 10px; display: flex; align-items: center; }
        .input-box { flex-grow: 1; padding: 12px; border: 1px solid #ccc; border-radius: 5px; font-size: 16px; }
        .send-button { padding: 12px 20px; background: #28a745; color: white; border: none; cursor: pointer; border-radius: 5px; margin-left: 10px; font-size: 16px; }
        .send-button:hover { background: #218838; }
        .bot-header { font-size: 24px; font-weight: bold; display: flex; align-items: center; justify-content: center; }
        .bot-header img { width: 30px; margin-right: 10px; }
        .message-container { display: flex; flex-direction: column; margin-bottom: 10px; }
        .user-label, .bot-label { font-weight: bold; margin-bottom: 2px; }
    </style>
</head>
<body>
    <div class="bot-header">
        <img src="https://cdn-icons-png.flaticon.com/512/4228/4228698.png" alt="Bot Icon"> 
        💊 Drug Discovery Chatbot
    </div>

    <div class="chat-container">
        <div class="chat-box" id="chat-box"></div>
        <div class="input-container">
            <input type="text" id="user-input" class="input-box" placeholder="Ask me anything about drug discovery..." />
            <button class="send-button" onclick="sendMessage()">Send</button>
        </div>
    </div>

    <script>
        async function sendMessage() {
            let userInput = document.getElementById("user-input").value.trim();
            if (!userInput) return;

            let chatBox = document.getElementById("chat-box");

            // Add user message with label "You"
            chatBox.innerHTML += "<div class='message-container'><div class='user-label'>You:</div><span class='user-message'>" + userInput + "</span></div>";
            document.getElementById("user-input").value = "";
            chatBox.scrollTop = chatBox.scrollHeight;

            let response = await fetch("/chat-ui/ask", {
                method: "POST",
                headers: { "Content-Type": "application/json" },
                body: JSON.stringify({ question: userInput })
            });

            let data = await response.json();
            let botResponse = formatResponse(data.answer);

            // Add bot message with label "Bot"
            chatBox.innerHTML += "<div class='message-container'><div class='bot-label'>MedXpert:</div><span class='bot-message'>" + botResponse + "</span></div>";
            chatBox.scrollTop = chatBox.scrollHeight;
        }

        function formatResponse(response) {
            return response.replace(/\\n/g, "<br>");  // Preserve line breaks in responses
        }
    </script>
</body>
</html>
"""

# Serve the HTML Chat UI
@app.get("/chat-ui", response_class=HTMLResponse)
async def chat_ui():
    return HTMLResponse(content=html_content)

# Chatbot API endpoint for the web UI
@app.post("/chat-ui/ask")
async def chat_with_gpt(request: Request):
    data = await request.json()
    user_question = data.get("question", "")

    if not user_question:
        return JSONResponse({"answer": "Please enter a valid question."})

    # Query OpenAI GPT
    try:
        response = openai.ChatCompletion.create(
            model="gpt-4",  # Change to "gpt-3.5-turbo" if needed
            messages=[
                {"role": "system", "content": "You are an AI assistant specializing in drug discovery."},
                {"role": "user", "content": user_question}
            ]
        )
        answer = response["choices"][0]["message"]["content"]

        # Format response (break paragraphs properly)
        formatted_answer = answer.replace("\n", "<br>")

    except Exception as e:
        formatted_answer = f"Error: {e}"

    return JSONResponse({"answer": formatted_answer})
