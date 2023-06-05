import json
def save(mapping, f):
    if (type(mapping) == "list"):
        with open(f, "w") as fi:
            f.write(*mapping, sep='\n')
    else:
        with open(f, "w") as fi:
            json.dump(mapping, fi, indent=2)  