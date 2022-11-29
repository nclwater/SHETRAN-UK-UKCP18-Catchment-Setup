import os
import pandas as pd
import shutil

root = "I:/SHETRAN_GB_2021/"
Copy_From = root + "05_Climate_Change_Simulations/UKCP18rcm_220708_APM_UK/"
Copy_To = root + "05_Climate_Change_Simulations/UKCP18rcm_Autocal_UDM_Baseline/"

RCPs = ['01', '04', "05", '06', '07', '08', '09', '10', '11', '12', '13', '15']

climate_files = ["_Cells.asc", "_PET.csv", "_Temp.csv", "_Precip.csv"]

# Make a list of catchments for copying if needed and read them in via a CSV file:
simulations_for_copying = pd.read_csv(Copy_From + "bc_01/UKCP18_runtime_estimates_bc_01.csv")
simulations_for_copying = simulations_for_copying["gauge_id"].astype(str)

# CEH_simulations = os.listdir(Copy_To + "bc_01")

# print(simulations_for_copying)
for rcp in RCPs:

    print(rcp)
    if not os.path.exists(Copy_To + "bc_" + rcp):
        os.mkdir(Copy_To + "bc_" + rcp)

    # Run through each folder:
    for sim in simulations_for_copying:

        simulation_folder = Copy_To + "bc_" + rcp + "/" + sim + "/"

        if not os.path.exists(simulation_folder):
            os.mkdir(simulation_folder)

        # Copy each of the climate files across (overwriting existing files):
        for cf in climate_files:

            src = Copy_From + "bc_" + rcp + "/" + sim + "/" + sim + cf
            dst = simulation_folder + sim + cf
            shutil.copy2(src, dst)
