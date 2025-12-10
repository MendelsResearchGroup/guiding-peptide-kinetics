import numpy as np

# YYDPETGTWY
long_to_short = {
    "chignolin":  "WT",
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
    "YYDPETGQWY": "T7Q", # Not enough grid size
    "YYDPETGRWY": "T7R",
    "YYDPETGYWY": "T7Y",
    "YYDPETGGWY": "T7G",
    "YYDPETGDWY": "T7D",
    "RYDPETGTWY": "Y0R",
    "QYDPETGTWY": "Y0Q",
    "EYDPETGTWY": "Y0E",
    "AYDPETGTWY": "Y0A",
    "YYDMETGTWY": "P3M",
    "YYDDETGTWY": "P3D",
    "YYDRETGTWY": "P3R",
    "YYDCETGTWY": "P3C",
    "YYDPERGTWY": "T5R",
    "YYDPEGGTWY": "T5G",
    "YYDPEYGTWY": "T5Y", # Not yet sampled
    "YYDPEDGTWY": "T5D",# Not yet sampled
    "YYDPATGTWY": "E4A",
    "YYDPGTGTWY": "E4G",
    "YYDPRTGTWY": "E4R",
    "YYDPYTGTWY": "E4Y", # Starts from high rmsd
    "YYKPETGTWY": "D2K",
    "YYYPETGTWY": "D2Y",
    # "YYDPETGTWK": "Y9K", outlier?
    
}

short_to_long = {v: k for k, v in long_to_short.items()}

short_to_medium = {
    "D2A": "Asp2Ala",
    "D2R": "Asp2Arg",
    "D2N": "Asp2Asn",
    "D2C": "Asp2Cys",
    "D2E": "Asp2Glu",
    "D2K": "Asp2Lys",
    "D2M": "Asp2Met",
    "D2Y": "Asp2Tyr",
    "P3R": "Pro3Arg",
    "P3D": "Pro3Asp",
    "P3C": "Pro3Cys",
    "P3M": "Pro3Met",
    "T7R": "Thr7Arg",
    "T7D": "Thr7Asp",
    "T7Q": "Thr7Gln",
    # "T7G": "Thr7Gly",
    "T7Y": "Thr7Tyr",
    "T7V": "Thr7Val",
    "Y0A": "Tyr0Ala",
    "Y0R": "Tyr0Arg",
    "Y0Q": "Tyr0Gln",
    "Y0E": "Tyr0Glu",
    "Y9A": "Tyr9Ala",
    "Y9R": "Tyr9Arg",
    "Y9Q": "Tyr9Gln",
    "Y9E": "Tyr9Glu",
    "Y9G": "Tyr9Gly",
    "Y9K": "Tyr9Lys",
    "Y9V": "Tyr9Val",
    "E4A": "Glu4Ala", 
    "E4G": "Glu4Gly", 
    "E4R": "Glu4Arg", 
    "E4Y": "Glu4Tyr"
}

groupByResidue = {
    0: ["Y0A", "Y0R", "Y0Q", "Y0E"],
    2: ["D2A", "D2C", "D2M", "D2N", "D2R", "D2E", "D2K", "D2Y"],
    3: ["P3C", "P3M", "P3R", "P3D"],
    4: ["E4A", "E4G", "E4R", "E4Y"],
    7: ["T7V", "T7Q", "T7R", "T7Y", 'T7G', 'T7D'],
    9: ["Y9E", "Y9G", "Y9Q", "Y9R", "Y9V", "Y9A"],
}

groupByProperty = {
    "RHKDE": ["D2R", "D2E", "Y9E", "Y9R", "T7R"],
    "STNQ": ["D2N", "Y9Q", "T7Q"],
    "CUGP": ["D2C", "Y9G", "T7G", "T7D"],
    "AVILMFYW": ["D2A", "D2M", "Y9A", "Y9V", "T7V", "T7Y", "Y0A"],
}

proteins = list(long_to_short.keys())

thresholds = np.arange(0.28, 0.5 + 1e-9, 0.02)
res_colors = {
    0: "red",
    2: "orange",
    3: "purple",
    4: "brown",
    5: "black",
    7: "blue",
    9: "green"
}
