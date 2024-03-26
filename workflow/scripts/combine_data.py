"""
This script will eventually do something.
"""

import os
import pandas as pd

# the working directory should house the analysis/ folder
working_dir = os.getcwd()

# create a dictionary that will contain all the data categories and their associated files
categories = {
    "don_hbonds_nr": [],
    "prot_don_hbonds": [],
    "acc_hbonds_nr": [],
    "deprot_acc_hbonds": [],
    "single_don_exo_amines": [],
    "dual_don_exo_amines": [],
    "don_hbonds_geom": [],
    "don_hbonds_geom_filtered": [],
    "nuc_data": []
}

# prepare a list of files from each category folder to combine
for cat in categories.keys():
    for file in os.listdir(working_dir + "/analysis/" + cat):
        # select files with the correct suffix
        if file[-(len(cat) + 5):] == "_" + cat + ".csv":
            # find out which files actually contain data
            contains_data = False
            with open(working_dir + "/analysis/" + cat + "/" + file, "r") as read_file:
                num_lines = 0
                for lines in read_file:
                    if num_lines > 0:
                        contains_data = True
                        break
                    else:
                        num_lines += 1
            # only include files that contain data
            if contains_data:
                categories[cat].append(working_dir + "/analysis/" + cat + "/" + file)

# create the combined/ folder to write csv files to
os.mkdir(working_dir + "/combined")

# create the combined dataframes and write to csv files
pd.concat([pd.read_csv(file) for file in categories["don_hbonds_nr"]]).to_csv("combined/don_hbonds_nr_c.csv", index=False)
pd.concat([pd.read_csv(file) for file in categories["prot_don_hbonds"]]).to_csv("combined/prot_don_hbonds_c.csv", index=False)
pd.concat([pd.read_csv(file) for file in categories["acc_hbonds_nr"]]).to_csv("combined/acc_hbonds_nr_c.csv", index=False)
pd.concat([pd.read_csv(file) for file in categories["deprot_acc_hbonds"]]).to_csv("combined/deprot_acc_hbonds_c.csv", index=False)
pd.concat([pd.read_csv(file) for file in categories["single_don_exo_amines"]]).to_csv("combined/single_don_exo_amines_c.csv", index=False)
pd.concat([pd.read_csv(file) for file in categories["dual_don_exo_amines"]]).to_csv("combined/dual_don_exo_amines_c.csv", index=False)
pd.concat([pd.read_csv(file) for file in categories["don_hbonds_geom"]]).to_csv("combined/don_hbonds_geom_c.csv", index=False)
pd.concat([pd.read_csv(file) for file in categories["don_hbonds_geom_filtered"]]).to_csv("combined/don_hbonds_geom_filtered_c.csv", index=False)
pd.concat([pd.read_csv(file) for file in categories["nuc_data"]]).to_csv("combined/nuc_data_c.csv", index=False)
