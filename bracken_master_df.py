#this script generates a master file than can be used for the analysis of the resulting bracken files
#it needs to be edited for the location of the data in input_dir and will output there under the name master_bracken_tidy.csv
#this script HAS to match the generated bracken files e.g the file nomenclature leads to its accurate classification under the spreadsheet
#             TO RUN:
#             python ./bracken_master_df.py

import os
import re
import pandas as pd

# Define folder containing all the Bracken output files, these can be exctracted frok the folders into a bracken specific folder
input_dir = "/home/s1887840/scratch/merge/controls/kraken"

# Output file name and path, can be changed if there is a need for a specific location
output_file = os.path.join(input_dir, "master_bracken_tidy.csv")

# Regex to extract metadata from filename, e.g. brkn_noplant_c05_622320_03.txt
# The pattern is flexible (e.g. bracken/brkn but not TOO flexible)
pattern = re.compile(r"^(bracken|brkn)_(\w+)_c(\d+)_([0-9]+)_([0-9]+)\.txt$")

#list to store dataframes
all_dfs = []

#loop over all .txt files in the directoy
for filename in os.listdir(input_dir):
    if filename.endswith(".txt"):
        match = pattern.match(filename)
        if not match:
            print(f"Skipping file with unexpected name format: {filename}")
            continue

        tool, database, confidence, sample_id, wave = match.groups()
        #tool, database, confidence, batch, barcode = match.groups()
        #for the batch analysis

        file_path = os.path.join(input_dir, filename)
        df = pd.read_csv(file_path, sep="\t")

        #Add metadata columns
        df["tool"] = tool
        df["database"] = database
        df["confidence"] = confidence
        df["sample_id"] = sample_id
        df["wave"] = wave
        #df["batch"] = batch
        #df["barcode"] = barcode
        #for the batch analysis

        all_dfs.append(df)

#combine all dataframes
if all_dfs:
    master_df = pd.concat(all_dfs, ignore_index=True)
    master_df.to_csv(output_file, index=False)
    print(f"\n✅ Master tidy file created: {output_file}")
    print(f"Total rows: {len(master_df)}")
# in case the folder has no files or something goes wrong
else:
    print("No valid files found in the directory.")
