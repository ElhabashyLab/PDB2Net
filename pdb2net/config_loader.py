import json
import os

CONFIG_PATH = os.path.join(os.path.dirname(__file__), "config.json")
_config_cache = None  # Cache für die geladene Konfiguration


def load_config():
    """
    Lädt die Konfigurationsdatei config.json mit Lazy Loading und Fehlerhandling.
    Falls die Datei nicht existiert oder fehlerhaft ist, wird eine Standardkonfiguration zurückgegeben.
    """
    global _config_cache
    if _config_cache is not None:
        return _config_cache  # Falls bereits geladen, verwende den Cache

    try:
        with open(CONFIG_PATH, "r") as f:
            _config_cache = json.load(f)  # Datei lesen und parsen
    except FileNotFoundError:
        print(f"⚠️ Warnung: {CONFIG_PATH} nicht gefunden! Es wird eine Standardkonfiguration verwendet.")
        _config_cache = {}
    except json.JSONDecodeError as e:
        print(f"❌ Fehler: {CONFIG_PATH} ist fehlerhaft! Bitte überprüfe die JSON-Syntax. Fehler: {e}")
        exit(1)

    return _config_cache


# Lazy Loading: Die Konfiguration wird erst bei Zugriff geladen
config = load_config()
