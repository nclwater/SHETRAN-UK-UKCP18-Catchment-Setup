# -----------------------------------------------------------------------------
# Setup UKCP18 SHETRAN Simulations
# 11/07/2022 & 04/10/2022
# Ben Smith (adapted from David Prichard)
# -----------------------------------------------------------------------------

"""
The following script has been taken from:
I:/SHETRAN_GB_2021/scripts/ukcp18rcm_210927/setup.py

It has been adapted so that instead of setting up separate control and future simulations, it instead uses a single
SHETRAN model to run through all the different periods. This is because we want to look at the periods in terms of
warming levels, rather than based on decades. As such, there are no "control" and "future" folders.

This also uses the APM soil/subsurface dataset, whereas the previous used the more basic setup.

This script will set up un bias controlled simulation by taking the historical simulations and adding climate data
from the UKCP18 datasets. Once you've done that then you need to run the bias_correction_quantile_mapping.py script.
That will make another setup, which can then be used for the Climate Change OpenCLIM analysis.

Once you have set up these catchments, check that the library file has the correct end date. If you set up the
historical catchments with different start/end dates to the ones referenced in this script then it will not correctly
match the toReplace/ReplacementStrings strings and it will not edit the simulation time and it will probably run for
30 years instead of 100.

"""

import os
import shutil
import numpy as np
import xarray as xr
import itertools
# import pandas as pd
# import sys
# import math
# from multiprocessing import Pool, Manager
# from time import sleep

# -----------------------------------------------------------------------------
# --- Setup User Inputs -------------------------------------------------------
# -----------------------------------------------------------------------------

root = 'I:/'
baseline_folder = root + 'SHETRAN_GB_2021/historical_220601_GB_APM/'
ukcp18_folder = root + 'UKCP18_UK_12km/'
ukcp18_pet_folder = root + 'SHETRAN_GB_2021/ukcp18rcm_pet/'
setup_folder = root + 'SHETRAN_GB_2021/UKCP18rcm_220708_APM_GB_NI/'
simulation_list = "I:/SHETRAN_GB_2021/scripts/UKCP18rcm_220708_APM_GB/Simulation_Setup_List.csv"
# models = ['01', '04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '15']
models = ['05']

make_folder_structure = True
periods = ['19801201-19901130', '19901201-20001130', '20001201-20101130', '20101201-20201130', '20201201-20301130',
           '20301201-20401130', '20401201-20501130', '20501201-20601130', '20601201-20701130', '20701201-20801130']


# -----------------------------------------------------------------------------
# --- Setup Functions ---------------------------------------------------------
# -----------------------------------------------------------------------------

def read_climate_data(root_folder, rcp, variable_name, period_list):

    first_loop = True

    # Run through the different decades, bolting the required catchment data into a common dataframe.
    for period in period_list:

        subpath = root_folder + "/" + variable_name + "/day/latest/" + variable_name + \
                  "_rcp85_land-rcm_uk_12km_" + rcp + "_day_" + period + ".nc"

        DS = xr.open_dataset(subpath)

        if first_loop:
            DS_all_periods = DS
            first_loop = False
        else:
            DS_all_periods = xr.merge([DS_all_periods, DS])

    return DS_all_periods


def get_variable(climate_xarray, variable, llx, lly, cellsize, ncols, nrows, outFileTS):
    """
    This will collect the different climate datasets.
    :param variable:
    :param llx:
    :param lly:
    :param cellsize:
    :param ncols:
    :param nrows:
    :param m: which of the Climate Models to use (1,4:15)
    :param scenario: A list of decadal dates for which there is model data.
    :param outFileTS:
    :return:
    """

    urx = llx + ncols * cellsize
    ury = lly + nrows * cellsize

    llx12 = int(llx / 12000) * 12000
    lly12 = int(lly / 12000) * 12000

    urx12 = (int(urx / 12000) + 1) * 12000
    ury12 = (int(ury / 12000) + 1) * 12000

    ds_subset = climate_xarray.sel(projection_y_coordinate=slice(lly12, ury12),
                                   projection_x_coordinate=slice(llx12, urx12))
    df = ds_subset[variable].to_dataframe()

    df = df.unstack(level=['projection_y_coordinate', 'projection_x_coordinate'])

    yCoords = list(df.columns.levels[1])
    yCoords.sort(reverse=True)

    xCoords = list(df.columns.levels[2])
    xCoords.sort(reverse=False)

    orderedDf = df.loc[:, list(itertools.product([variable], yCoords, xCoords))]

    orderedDf.to_csv(outFileTS, index=False, header=np.arange(1, len(orderedDf.columns) + 1))


def setItUp(catchment_name, m, rainfall, temperature, pet, q=None):

    setup_subfolder = setup_folder + m + "/" + catchment_name + "/"

    files_to_copy = [
        catchment_name + "_DEM.asc",
        catchment_name + "_Lake.asc",
        catchment_name + "_LandCover.asc",
        catchment_name + "_LibraryFile.xml",
        catchment_name + "_Mask.asc",
        catchment_name + "_MinDEM.asc",
        catchment_name + "_Soil.asc"
    ]
    for file in files_to_copy:
        shutil.copy(baseline_folder + catchment_name + "/" + file, setup_subfolder + file)

    with open(setup_subfolder + catchment_name + "_Mask.asc", "r") as mask_file:
        ncols = int(mask_file.readline().rstrip().split()[1])
        nrows = int(mask_file.readline().rstrip().split()[1])
        llx = float(mask_file.readline().rstrip().split()[1])
        lly = float(mask_file.readline().rstrip().split()[1])
        cellsize = float(mask_file.readline().rstrip().split()[1])

    # ---

    # Read in the library file from the APM simulation.
    with open(setup_subfolder + catchment_name + "_LibraryFile.xml", 'r') as file:
        library_file = file.read()

    # Edit the library file to reflect the changed time period:
    toReplace = "<StartDay>1</StartDay>\n<StartMonth>1</StartMonth>\n<StartYear>1980</StartYear>\n<EndDay>1</EndDay" \
                ">\n<EndMonth>1</EndMonth>\n<EndYear>2011</EndYear>"

    replacementText = "<StartDay>01</StartDay>\n<StartMonth>12</StartMonth>\n<StartYear>1980</StartYear>\n<EndDay>30" \
                      "</EndDay>\n<EndMonth>11</EndMonth>\n<EndYear>2080</EndYear> "

    # Replace the target string:
    library_file = library_file.replace(toReplace, replacementText)

    # Write the library file out again
    with open(setup_subfolder + catchment_name + "_LibraryFile.xml", 'w') as file:
        file.write(library_file)

    # ---

    # Make climate input time series files
    print("----- Writing Rainfall")
    get_variable(rainfall, 'pr', llx, lly, cellsize, ncols, nrows,
                 setup_subfolder + "/" + catchment_name + "_Precip.csv")

    print("----- Writing Temperature")
    get_variable(temperature, 'tas', llx, lly, cellsize, ncols, nrows,
                 setup_subfolder + "/" + catchment_name + "_Temp.csv")

    print("----- Writing PET")
    get_variable(pet, 'pet', llx, lly, cellsize, ncols, nrows,
                 setup_subfolder + catchment_name + "_PET.csv")

    # ---

    print("----- Completing")
    idDict = {}
    ticker = 1
    newMap = np.arange(1, ncols * nrows + 1).reshape((nrows, ncols))

    for j in range(nrows):
        for i in range(ncols):
            xc = int(((i * cellsize) + llx) / 12000)

            yc = int((((nrows - 1 - j) * cellsize) + lly) / 12000)

            id = str(xc) + "," + str(yc)

            if id not in idDict.keys():
                idDict[id] = ticker
                ticker += 1
            else:
                pass

            newMap[j][i] = idDict[id]

    head = 'ncols\t' + str(ncols) + '\nnrows\t' + str(nrows) + '\nxllcorner\t' + str(llx) + '\nyllcorner\t' + str(
        lly) + '\ncellsize\t' + str(cellsize) + '\nNODATA_value\t-9999'

    np.savetxt(setup_subfolder + catchment_name + "_Cells.asc", newMap,
               delimiter=' ', header=head, fmt='%.0f', comments='')


# -----------------------------------------------------------------------------
# --- Create a list of catchments to process ----------------------------------
# -----------------------------------------------------------------------------

catchments = []
with open(simulation_list, 'r') as fh:
    fh.readline()
    for line in fh:
        catchments.append(line.rstrip().split(',')[0])
# catchments = ["201002", "201006"]

# -----------------------------------------------------------------------------
# --- Create the necessary folder structure if needed -------------------------
# -----------------------------------------------------------------------------

if make_folder_structure:
    RCP_root_folder = setup_folder + "/"
    if not os.path.exists(RCP_root_folder):
        print("Making RCP root folder...")
        os.mkdir(RCP_root_folder)

    for m in models:
        model_path = RCP_root_folder + m + "/"
        if not os.path.exists(model_path):
            print("Making RCP Simulation folder...")
            os.mkdir(model_path)
        for catch in catchments:
            catch_path = model_path + catch + '/'
            if not os.path.exists(catch_path):
                print("Making simulation folder...")
                os.makedirs(catch_path)


# -----------------------------------------------------------------------------
# --- Setup the catchment files -----------------------------------------------
# -----------------------------------------------------------------------------

if __name__ == '__main__':

    for rcp in models:
        print(rcp)

        print("  Reading rainfall...")
        rainfall_dataset = read_climate_data(root_folder=ukcp18_folder + "/" + rcp,
                                             rcp=rcp, variable_name='pr',
                                             period_list=periods)

        print("  Reading temperature...")
        temperature_dataset = read_climate_data(root_folder=ukcp18_folder + "/" + rcp,
                                                rcp=rcp, variable_name='tas',
                                                period_list=periods)

        print("  Reading PET...")
        pet_dataset = read_climate_data(root_folder=ukcp18_pet_folder + "/" + rcp,
                                        rcp=rcp, variable_name='pet',
                                        period_list=periods)

        for catch in catchments:
            print("---", catch)
            setItUp(catchment_name=catch, m=rcp,
                    rainfall=rainfall_dataset, temperature=temperature_dataset, pet=pet_dataset)
