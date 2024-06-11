import csv
import os

from biothings import config

logging = config.logger


def load_interactions(data_folder: str):
    infile = os.path.join(data_folder, "approved_drug_detailed_interactions.csv")
    assert os.path.exists(infile)

    # pk: Ligand ID,Target ID,Target Ligand ID,Target Species
    # inner join of primary_targets_csv[pk] and detailed_csv[pk] is primary_targets_csv[pk]
    with open(infile, "r") as data_fd:
        next(data_fd)  # skip metadata row
        reader = csv.DictReader(data_fd)
        idx = 0

        for row in reader:
            row["_id"] = str(idx)
            idx += 1
            yield row
