"""
This script will eventually do something.
"""

import os
import pandas as pd

# The working directory should house the analysis/ folder.
working_dir = os.getcwd()

# Create a dictionary that will contain all the data categories and their associated files.
categories = {
    "nuc_data": [],
    "b_factor_data": [],
    "don_hbonds_nr": [],
    "prot_don_hbonds_nr": [],
    "acc_hbonds_nr": [],
    "hbond_data": []
}

# Prepare a list of files from each category folder to combine.
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

# Create the combined/ folder to write csv files to.
os.mkdir(working_dir + "/combined")

# Create the combined dataframes and write to csv files.
pd.concat([pd.read_csv(file, na_filter=False) for file in categories["nuc_data"]]).to_csv("combined/nuc_data_c.csv", index=False)
pd.concat([pd.read_csv(file, na_filter=False) for file in categories["b_factor_data"]]).to_csv("combined/b_factor_data_c.csv", index=False)
pd.concat([pd.read_csv(file, na_filter=False) for file in categories["don_hbonds_nr"]]).to_csv("combined/don_hbonds_nr_c.csv", index=False)
pd.concat([pd.read_csv(file, na_filter=False) for file in categories["prot_don_hbonds_nr"]]).to_csv("combined/prot_don_hbonds_nr_c.csv", index=False)
pd.concat([pd.read_csv(file, na_filter=False) for file in categories["acc_hbonds_nr"]]).to_csv("combined/acc_hbonds_nr_c.csv", index=False)
pd.concat([pd.read_csv(file, na_filter=False) for file in categories["hbond_data"]]).to_csv("combined/hbond_data_c.csv", index=False)
