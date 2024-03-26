"""
This module will eventually do something.
"""

import sys
import os
import pandas as pd

# the working directory should house the analysis/ folder
working_dir = os.getcwd()

# set the directory pointing to the don_hbonds_geom/ folder
don_hbonds_geom_dir = working_dir + "/analysis/don_hbonds_geom"

# prepare a list of files from the don_hbonds_geom/ folder to combine
don_hbonds_geom_files = []
for file in os.listdir(don_hbonds_geom_dir):
    # select files with the correct suffix
    if file[-20:] == "_don_hbonds_geom.csv":
        # find out which files actually contain data
        contains_data = False
        with open(don_hbonds_geom_dir + "/" + file, "r") as read_file:
            num_lines = 0
            for lines in read_file:
                if num_lines > 0:
                    contains_data = True
                    break
                else:
                    num_lines += 1
        # only include files that contain data
        if contains_data:
            don_hbonds_geom_files.append(don_hbonds_geom_dir + "/" + file)

# prepare the new dataframe and write it to a csv file
don_hbonds_geom = pd.concat([pd.read_csv(file) for file in don_hbonds_geom_files])
don_hbonds_geom.to_csv("don_hbonds_geom_combined.csv", index=False)

# set the directory pointing to the don_hbonds_geom_filtered/ folder
don_hbonds_geom_filtered_dir = working_dir + "/analysis/don_hbonds_geom_filtered"

# prepare a list of files from the don_hbonds_geom_filtered/ folder to combine
don_hbonds_geom_filtered_files = []
for file in os.listdir(don_hbonds_geom_filtered_dir):
    # select files with the correct suffix
    if file[-29:] == "_don_hbonds_geom_filtered.csv":
        # find out which files actually contain data
        contains_data = False
        with open(don_hbonds_geom_filtered_dir + "/" + file, "r") as read_file:
            num_lines = 0
            for lines in read_file:
                if num_lines > 0:
                    contains_data = True
                    break
                else:
                    num_lines += 1
        # only include files that contain data
        if contains_data:
            don_hbonds_geom_filtered_files.append(don_hbonds_geom_filtered_dir + "/" + file)

# prepare the new dataframe and write it to a csv file
don_hbonds_geom_filtered = pd.concat([pd.read_csv(file) for file in don_hbonds_geom_filtered_files])
don_hbonds_geom_filtered.to_csv("don_hbonds_geom_filtered_combined.csv", index=False)
