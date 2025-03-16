import json
import os

# Define the path to the configuration file
CONFIG_PATH = os.path.join(os.path.dirname(__file__), "config.json")

# Global cache for configuration
_config_cache = None


def load_config():
    """
    Loads the configuration file (config.json) using lazy loading and error handling.

    If the file does not exist or contains errors, a default empty configuration is returned.

    Returns:
        dict: The loaded configuration dictionary.
    """
    global _config_cache
    if _config_cache is not None:
        return _config_cache  # Return cached configuration if already loaded

    try:
        with open(CONFIG_PATH, "r") as f:
            _config_cache = json.load(f)  # Read and parse the JSON file
    except FileNotFoundError:
        print(f"Warning: {CONFIG_PATH} not found! Using default configuration.")
        _config_cache = {}
    except json.JSONDecodeError as e:
        print(f"Error: {CONFIG_PATH} contains invalid JSON syntax. Please check the file. Details: {e}")
        exit(1)  # Terminate the program in case of a JSON syntax error

    return _config_cache


# Load the configuration at module import (lazy loading)
config = load_config()
