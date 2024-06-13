import csv
import os

from biothings import config

logging = config.logger


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
            yield row


def load_ligands(data_folder: str):
    ligands_file = os.path.join(data_folder, "ligands.csv")
    assert os.path.exists(ligands_file)

    with open(ligands_file, "r") as data_fd:
        next(data_fd)  # skip metadata row
        reader = csv.DictReader(data_fd)

        for row in reader:
            row["_id"] = row["Ligand ID"]
            yield row
