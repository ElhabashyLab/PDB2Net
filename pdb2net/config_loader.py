import json
import os

# Load configuration from config.json
CONFIG_PATH = os.path.join(os.path.dirname(__file__), "config.json")

def load_config():
    """Loads the configuration from config.json"""
    with open(CONFIG_PATH, "r") as f:
        return json.load(f)

config = load_config()
