# ----------------------------------------------------------------
# -- 3bi Recreate SHETRAN Simulations for Extracting Soil Moisture
# ----------------------------------------------------------------
# -- Ben Smith ---------------------------------------------------
# -- 18/10/2022 --------------------------------------------------
# ----------------------------------------------------------------
# This script takes the existing UKCP18 simulations and copies their input files to a new repository.
# The library files are then altered to reduce the simulation duration to only around 30 years.
# This is to allow the simulations to be rerun to create soil moisture values that can be used by the CEH team.
# Previous simulations only produced soil moisture for the top 6 cells (as this was 2m in the autocalibration
# setup). In this setup, which uses the default SHETRAN prepare script, the subsurface layers are thinner, and
# so more layers are required to calculate the soil moisture of the top 2m (as requested by CEH).

# This is quite slow - I think due to the loading / editing / writing the csv's - but there.

# --- PREAMBLE

import os
import pandas as pd
import shutil

root = "I:/SHETRAN_GB_2021/05_Climate_Change_Simulations/UKCP18rcm_220708_APM_UK/"
root_CEH = "I:/SHETRAN_GB_2021/05_Climate_Change_Simulations/UKCP18rcm_181022_APM_UK_CEH/"

# models = ['01', '04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '15']
models = ['05']

# List out the catchment names:
catchments = pd.read_csv(root_CEH + "Catchment Ranking for Soil Moisture.csv")
catchments = catchments[catchments['Rank'] > 0]
catchments = catchments["Catchment"].to_list()
catchments = [str(c) for c in catchments]

simulation_statics = ['Cells.asc', 'DEM.asc', 'Lake.asc', 'LandCover.asc', 'Mask.asc', 'MinDEM.asc', 'Soil.asc']
simulation_climate = ['PET.csv', 'Precip.csv', 'Temp.csv']


# --- PROCESS RCPS:
for m in models:

    print(m)

    # List the subfolders:
    master_folder = root + "bc_" + m + "/"
    folders = os.listdir(master_folder)
    folders = [f for f in folders if "." not in f]
    folders = [f for f in folders if f in catchments]

    # Run through subfolders (i.e. catchments) and copy simulation inputs across:
    if not os.path.isdir(root_CEH + "bc_" + m):
        os.mkdir(root_CEH + "bc_" + m)

    for folder in folders:

        # Check whether the catchment has been done already:
        if os.path.isfile(root_CEH + "/bc_" + m + "/" + folder + "/" + folder + "_LibraryFile.xml"):
            print(" -- ", folder, "- already made")
            continue

        print(" -- ", folder)

        # Create a catchment folder if needed:
        if not os.path.isdir(root_CEH + "bc_" + m + "/" + folder + "/"):
            os.mkdir(root_CEH + "bc_" + m + "/" + folder)

        try:

            for file in simulation_statics:
                shutil.copy2(root + "bc_" + m + "/" + folder + "/" + folder + "_" + file,
                             root_CEH + "bc_" + m + "/" + folder + "/" + folder + "_" + file)

            # Do the same for climate files, but crop them down to 31 years (just to be safe) in the process:
            for file in simulation_climate:
                climate = pd.read_csv(root + "bc_" + m + "/" + folder + "/" + folder + "_" + file, nrows=360*31)
                climate.to_csv(root_CEH + "bc_" + m + "/" + folder + "/" + folder + "_" + file, index=False)

            # Now edit the library files:
            library_path = m + "/" + folder + "/" + folder + "_LibraryFile.xml"

            # Read in the library files:
            with open(root + "/bc_" + library_path, 'r') as file:
                library_file = file.read()

            # Edit the library file to reflect the changed time period:
            library_file = library_file.replace("<EndYear>2080</EndYear>", "<EndYear>2010</EndYear>")

            # Write the library file out again
            with open(root_CEH + "bc_" + library_path, 'w') as file:
                file.write(library_file)

        except Exception as e:
            print(folder.upper(), e)
