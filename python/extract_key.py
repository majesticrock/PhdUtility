import gzip

def extract_key(file_path, key):
    value = None
    if file_path.endswith(".gz"):
        with gzip.open(file_path, "rt") as file:
            for line in file:
                if line.startswith(f"# {key}:"):
                    value = float(line.split(":")[1].strip())
                    break
    else:
        with open(file_path, 'r') as file:
            for line in file:
                if line.startswith(f"# {key}:"):
                    value = float(line.split(":")[1].strip())
                    break
    return value