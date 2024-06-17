import os

import numpy as np
import pandas as pd
from biothings import config
from biothings.utils.dataload import dict_convert, dict_sweep

logging = config.logger

VAL_MAP = {"yes": True, "no": False}
process_key = lambda key: key.replace(" ", "_").lower()
process_val = lambda val: VAL_MAP[val] if type(val) == str and val in VAL_MAP.keys() else val


def preprocess_data(d: dict):
    d = dict_sweep(d, vals=["", np.nan], remove_invalid_list=True)
    d = dict_convert(d, keyfn=process_key)
    d = dict_convert(d, valuefn=process_val)
    return d


def load_ligands(data_folder: str):
    # pk: Ligand ID,Target ID,Target Ligand ID,Target Species
    # inner join of primary_targets_csv[pk] and detailed_csv[pk] is primary_targets_csv[pk]
    interactions_file = os.path.join(data_folder, "approved_drug_detailed_interactions.csv")
    ligands_file = os.path.join(data_folder, "ligands.csv")
    assert os.path.exists(interactions_file) and os.path.exists(ligands_file)

    ligands_df = (
        pd.read_csv(ligands_file, skiprows=1, dtype=object)
        .rename({"Ligand ID": "_id"}, axis=1)
        .set_index("_id")
    )
    interactions_df = pd.read_csv(interactions_file, skiprows=1, dtype=object).rename(
        {"Ligand": "Name"}, axis=1
    )

    # drop redundant common columns from interactions
    cols_to_drop = ligands_df.columns.intersection(interactions_df.columns)
    interactions = interactions_df.drop(cols_to_drop, axis=1).to_dict(orient="records")
    ligands = ligands_df.to_dict(orient="index")
    assert type(list(ligands.keys())[0]) == str

    for row in interactions:
        ligand_id = str(row["Ligand ID"])
        row.pop("Ligand ID")  # redundant since present in ligands

        # NOTE: we assume ligand IDs in interactions will be found in ligands
        if "interactions" not in ligands[ligand_id].keys():
            ligands[ligand_id]["interactions"] = []
        ligands[ligand_id]["interactions"].append(preprocess_data(row))

    for k, ligand in ligands.items():
        ligand["_id"] = k
        yield preprocess_data(ligand)
