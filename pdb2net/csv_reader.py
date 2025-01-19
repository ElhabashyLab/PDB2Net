import csv

def read_csv(file_path):
    """
    Liest eine CSV-Datei und extrahiert die Dateipfade.

    :param file_path: Pfad zur CSV-Datei
    :return: Liste von Dateipfaden
    """
    paths = []
    with open(file_path, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        if "file_path" not in reader.fieldnames:
            raise ValueError("Die CSV-Datei muss eine Spalte 'file_path' enthalten.")
        for row in reader:
            paths.append(row["file_path"])
    return paths
