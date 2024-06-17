import csv
import os

from biothings import config
from biothings.utils.dataload import dict_convert, dict_sweep

logging = config.logger

VAL_MAP = {"yes": True, "no": False}
process_key = lambda key: key.replace(" ", "_").lower()
process_val = lambda val: VAL_MAP[val] if val in VAL_MAP.keys() else val


def load_interactions(data_folder: str):
    # pk: Ligand ID,Target ID,Target Ligand ID,Target Species
    # inner join of primary_targets_csv[pk] and detailed_csv[pk] is primary_targets_csv[pk]
    interactions_file = os.path.join(data_folder, "approved_drug_detailed_interactions.csv")
    assert os.path.exists(interactions_file)

    with open(interactions_file, "r") as data_fd:
        next(data_fd)  # skip metadata row
        reader = csv.DictReader(data_fd)

        for row in reader:
            row["_id"] = row["Ligand ID"]

            # remove keys with empty string OR list of empty strings
            row = dict_sweep(row, vals=[""], remove_invalid_list=True)
            row = dict_convert(row, keyfn=process_key)
            row = dict_convert(row, valuefn=process_val)
            yield row


def load_ligands(data_folder: str):
    ligands_file = os.path.join(data_folder, "ligands.csv")
    assert os.path.exists(ligands_file)

    with open(ligands_file, "r") as data_fd:
        next(data_fd)  # skip metadata row
        reader = csv.DictReader(data_fd)

        for row in reader:
            row["_id"] = row["Ligand ID"]

            # remove keys with empty string OR list of empty strings
            row = dict_sweep(row, vals=[""], remove_invalid_list=True)
            row = dict_convert(row, keyfn=process_key)
            row = dict_convert(row, valuefn=process_val)
            yield row
