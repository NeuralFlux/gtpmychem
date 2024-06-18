import re

import numpy as np
from biothings.utils.dataload import dict_convert, dict_sweep

VAL_MAP = {"yes": True, "no": False}
process_key = lambda key: key.replace(" ", "_").lower()
process_val = lambda val: VAL_MAP[val] if isinstance(val, str) and val in VAL_MAP.keys() else val
intrs_rename_dict = {
    "Target Ensembl Gene ID": "Ensembl Gene",
    "Target Entrez Gene ID": "Entrez Gene",
    "Target Gene Name": "Gene Name",
    "Target Species": "Species",
}


def preprocess_ligands(d: dict):
    if isinstance(d["Synonyms"], str):
        d["Synonyms"] = d["Synonyms"].split("|")
    d = dict_sweep(d, vals=["", np.nan], remove_invalid_list=True)
    d = dict_convert(d, keyfn=process_key)
    d = dict_convert(d, valuefn=process_val)
    return d


def preprocess_intrs(d: dict):
    d["Name"] = d["Target"]
    if isinstance(d["Species"], str):
        d["Species"] = d["Species"].lower()
    if isinstance(d["Name"], str):
        d["Name"] = re.sub(r"</?sub>", "", d["Name"])  # replace subscript tags

    # redundant since present in ligands
    cols_to_drop = [
        "Ligand ID",
        "CAS Number",
        "Clinical Use Comment",
        "Ligand Synonyms",
        "Target",
        "Ligand",
        "Type",
        "SMILES",
    ]
    for col in cols_to_drop:
        d.pop(col)

    d = dict_sweep(d, vals=["", np.nan], remove_invalid_list=True)
    d = dict_convert(d, keyfn=process_key)
    return d
