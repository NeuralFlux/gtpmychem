import os

import pandas as pd
from biothings import config

import utils as ligand_utils

logging = config.logger


def load_ligands(data_folder: str):
    # pk: Ligand ID,Target ID,Target Ligand ID,Target Species
    # inner join of primary_targets_csv[pk] and detailed_csv[pk] is primary_targets_csv[pk]
    interactions_file = os.path.join(data_folder, "approved_drug_detailed_interactions.csv")
    ligands_file = os.path.join(data_folder, "ligands.csv")
    assert os.path.exists(interactions_file) and os.path.exists(ligands_file)

    ligands = (
        pd.read_csv(ligands_file, skiprows=1, dtype=object)
        .rename({"Ligand ID": "_id"}, axis=1)
        .set_index("_id")
        .to_dict(orient="index")
    )
    interactions = (
        pd.read_csv(interactions_file, skiprows=1, dtype=object)
        .rename(ligand_utils.intrs_rename_dict, axis=1)
        .to_dict(orient="records")
    )
    assert type(list(ligands.keys())[0]) == str

    for row in interactions:
        ligand_id = str(row["Ligand ID"])

        # NOTE: we assume ligand IDs in interactions will be found in ligands
        if "interaction_targets" not in ligands[ligand_id].keys():
            ligands[ligand_id]["interaction_targets"] = []
            ligands[ligand_id]["CAS Number"] = row["CAS Number"]
            ligands[ligand_id]["Clinical Use Comment"] = row["Clinical Use Comment"]

        ligands[ligand_id]["interaction_targets"].append(ligand_utils.preprocess_intrs(row))

    for k, ligand in ligands.items():
        ligand["_id"] = k
        yield ligand_utils.preprocess_ligands(ligand)
