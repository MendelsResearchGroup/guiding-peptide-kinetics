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
    "chignolin": "Wildtype",
}


groupByResidue = {
    2: ["D2A", "D2C", "D2M", "D2N", "D2R", "D2E"],
    7: ["T7V", "T7Q", "T7R", "T7Y"],
    9: ["Y9E", "Y9G", "Y9Q", "Y9R", "Y9V", "Y9A"],
}

groupByProperty = {
    "RHKDE": ["D2R", "D2E", "Y9E", "Y9R", "T7R"],
    "STNQ": ["D2N", "Y9Q", "T7Q"],
    "CUGP": ["D2C", "Y9G"],
    "AVILMFYW": ["D2A", "D2M", "Y9A", "Y9V", "T7V", "T7Y"],
}

proteins = [
    "chignolin",
    "YYCPETGTWY",
    "YYDPETGTWE",
    "YYRPETGTWY",
    "YYAPETGTWY",
    "YYDPETGQWY",
    "YYDPETGTWG",
    "YYDPETGTWQ",
    "YYNPETGTWY",
    "YYEPETGTWY",
    "YYDPETGTWR",
    "YYDPETGVWY",
    "YYDPETGYWY",
    "YYDPETGTWV",
    "YYDPETGRWY",
]

thresholds = np.arange(0.1, 0.50, 0.04)
