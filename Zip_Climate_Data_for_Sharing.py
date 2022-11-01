# ------------------------------
# Zip Climate Data for Sharing
# Ben Smith
# 21/10/2022
# ------------------------------
# This notebook runs through the UKCP18 SHETRAN simulations and extracts and zips the climate data so that it can be
# shared with the HBV team. It will also note all instances where not all the climate data is present, so that it
# can be recreated and rerun.
# ------------------------------


# --- PREAMBLE ---
import os
from zipfile import ZipFile, ZIP_DEFLATED
import sys
import pandas as pd
# import tarfile

root = "I:/SHETRAN_GB_2021/"
CC = root + "05_Climate_Change_Simulations/UKCP18rcm_220708_APM_UK/"
zipfolder = CC + "/Bias Controlled Climate Zipfiles/"

# Second pass - path to additional files for zipping:
only_zip_some_folders = True
if only_zip_some_folders:
    files_for_zip = pd.read_csv(zipfolder + "Missing_Climate_Files_for_Copying.csv")

# models = ["01"]  # ["01", '04', "05", '06', '07', '08', '09', '10', '11', '12', '13', '15']
models = [sys.argv[1]]  # For running via command line, add the model name after the script name.

climate_files = ["_Cells.asc", "_PET.csv", "_Temp.csv", "_Precip.csv"]
# --- PROCESS ---

# Run through RCP models:
for m in models:

    # Create list to store details of catchments with missing data:
    catchment_issues = []

    # Create a folder to take the zipped files for that model:
    master_zip = zipfolder + "bc_" + m + '/'
    if not os.path.isdir(master_zip):
        os.mkdir(master_zip)

    # List the subfolders top run through:
    master_folder = CC + "bc_" + m + "/"

    # If you are only zipping some folders, take these from the list:
    if only_zip_some_folders:
        folders = files_for_zip.loc[files_for_zip['RCP'] == int(m)]["Catchment"]

    # Else list them from the available folders in the directory:
    else:
        folders = os.listdir(master_folder)
        # Remove any files from the list of folders:
        folders = [f for f in folders if "." not in f]

    # Run through subfolders (i.e. catchments):
    for folder in folders:

        print(m, "-", folder)

        # Set up a zipfile for holding the climate data for that catchment:
        with ZipFile(master_zip + "/" + folder + '.zip', 'w', ZIP_DEFLATED) as zip:
        # with tarfile.open(master_zip + "/" + folder + '.tar.gz', "w:gz") as tar:

            # Run through the four climate files:
            for c in climate_files:

                # make a path to the file to be copied:
                clim_file = master_folder + folder + "/" + folder + c

                # Check that the file is there and note if it is not:
                if not os.path.exists(clim_file):
                    catchment_issues.append(m + "_" + folder + '-' + c)
                    continue

                # Check that the file has data in it and note if it does not:
                if os.stat(clim_file).st_size < 1:
                    catchment_issues.append(m + "_" + folder + '-' + c)
                else:
                    # If the file exists and contains data, write it to a zip file:
                    zip.write(clim_file, folder + c)
                    # tar.add(clim_file)

    print(catchment_issues)
    with open(zipfolder + m + '_climate_data_summary_second_pass.txt', 'w') as f:
        for i in catchment_issues:
            f.write(i)
