import numpy as np

# YYDPETGTWY
mutation_map = {
    "YYAPETGTWY": "D2A",
    "YYCPETGTWY": "D2C",
    "YYMPETGTWY": "D2M",
    "YYNPETGTWY": "D2N",
    "YYRPETGTWY": "D2R",
    "YYEPETGTWY": "D2E",
    "YYDPETGTWE": "Y9E",
    "YYDPETGTWG": "Y9G",
    "YYDPETGTWQ": "Y9Q",
    "YYDPETGTWR": "Y9R",
    "YYDPETGTWV": "Y9V",
    "YYDPETGTWA": "Y9A",
    "YYDPETGVWY": "T7V",
    "YYDPETGQWY": "T7Q",
    "YYDPETGRWY": "T7R",
    "YYDPETGYWY": "T7Y",
    "YYDPETGGWY": "T7G",
    "YYDPETGDWY": "T7D",
    "chignolin":  "WT",
    "AYDPETGTWY": "Y0A",
    "YYDCETGTWY": "P3C",
    "RYDPETGTWY": "Y0R",
    "QYDPETGTWY": "Y0Q",
    "EYDPETGTWY": "Y0E",
    "YYDMETGTWY": "P3M",
    "YYDDETGTWY": "P3D",
    "YYDRETGTWY": "P3R"
}


groupByResidue = {
    0: ["Y0A", "Y0R", "Y0Q", "Y0E"],
    2: ["D2A", "D2C", "D2M", "D2N", "D2R", "D2E"],
    3: ["P3C", "P3M", "P3R", "P3D"],
    7: ["T7V", "T7Q", "T7R", "T7Y", 'T7G', 'T7D'],
    9: ["Y9E", "Y9G", "Y9Q", "Y9R", "Y9V", "Y9A"],
}

groupByProperty = {
    "RHKDE": ["D2R", "D2E", "Y9E", "Y9R", "T7R"],
    "STNQ": ["D2N", "Y9Q", "T7Q"],
    "CUGP": ["D2C", "Y9G", "T7G", "T7D"],
    "AVILMFYW": ["D2A", "D2M", "Y9A", "Y9V", "T7V", "T7Y", "Y0A"],
}

proteins = list(mutation_map.keys())

thresholds = np.arange(0.1, 0.50 + 1e-9, 0.04)
res_colors = {
    0: "red",
    2: "orange",
    3: "purple",
    7: "blue",
    9: "green"
}
