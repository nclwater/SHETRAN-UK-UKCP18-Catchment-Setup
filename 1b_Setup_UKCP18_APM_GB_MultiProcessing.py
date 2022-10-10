# -----------------------------------------------------------------------------
# Setup UKCP18 SHETRAN Simulations
# 11/07/2022
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

"""

import os
import shutil
import xarray as xr
import itertools
import pandas as pd
import numpy as np

# import math
from multiprocessing import Pool, Manager
# from time import sleep

# -----------------------------------------------------------------------------
# --- Setup User Inputs -------------------------------------------------------
# -----------------------------------------------------------------------------

root = 'I:/'
baseline_folder = root + 'SHETRAN_GB_2021/historical_220601_GB_APM/'
ukcp18_folder = root + 'UKCP18_UK_12km/'
ukcp18_pet_folder = root + 'SHETRAN_GB_2021/ukcp18rcm_pet/'
setup_folder = root + 'SHETRAN_GB_2021/UKCP18rcm_220708_APM_GB_NI/'
# simulation_list = root + 'SHETRAN_GB_2021/inputs/CAMELS_GB_topographic_attributes.csv'
# Simulation List for Northern Ireland additions:
simulation_list = "C:/Users/nbs65/OneDrive - Newcastle University/Python Code/SHETRAN_generic_catchment_setup/Simulation_Setup_List.csv"
models = ['01', '04'] #  , '05', '06', '07', '08', '09', '10', '11', '12', '13', '15']

make_folder_structure = True

n_processes = 10


# -----------------------------------------------------------------------------
# --- Setup Functions ---------------------------------------------------------
# -----------------------------------------------------------------------------

# This is our callback function. When a worker (i.e. a python process running runMe()) has FINISHED, this
# function will be called with whatever was returned (line 34, return taskNumber).
def workerDone(catch):
    # We print a success message here
    print(" -", catch, " completed")
    pass


def get_variable(variable, llx, lly, cellsize, ncols, nrows, m, scenario, outFileTS):
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

    # Run through the different decades, bolting the required catchment data into a common dataframe.
    orderedDfs = []
    for period in scenario:

        subpath = (
                m + "/" + variable + "/day/latest/" + variable
                + "_rcp85_land-rcm_uk_12km_" + m + "_day_" + period + ".nc"
        )

        if variable == 'pet':
            DS = xr.open_dataset(ukcp18_pet_folder + subpath)
        else:
            DS = xr.open_dataset(ukcp18_folder + subpath)

        ds_subset = DS.sel(projection_y_coordinate=slice(lly12, ury12),
                           projection_x_coordinate=slice(llx12, urx12))
        df = ds_subset[variable].to_dataframe()

        df = df.unstack(level=['projection_y_coordinate', 'projection_x_coordinate'])

        yCoords = list(df.columns.levels[1])
        yCoords.sort(reverse=True)

        xCoords = list(df.columns.levels[2])
        xCoords.sort(reverse=False)

        orderedDf = df.loc[:, list(itertools.product([variable], yCoords, xCoords))]
        orderedDfs.append(orderedDf)

    all_periods = pd.concat(orderedDfs)  # .sort_index()#.loc[startime:endtime]

    all_periods.to_csv(outFileTS, index=False, header=np.arange(1, len(all_periods.columns) + 1))


def setItUp(catch, m, q=None):  # controlOrFuture
    # print(m, controlOrFuture, catch)

    setup_subfolder = setup_folder + m + "/" + catch + "/"
    #  if not os.path.exists(newPath):
    #      os.makedirs(newPath)

    # print(setup_subfolder)

    files_to_copy = [
        catch + "_DEM.asc",
        catch + "_Lake.asc",
        catch + "_LandCover.asc",
        catch + "_LibraryFile.xml",
        catch + "_Mask.asc",
        catch + "_MinDEM.asc",
        catch + "_Soil.asc"
    ]
    for file in files_to_copy:
        # print(file)
        shutil.copy(baseline_folder + catch + "/" + file, setup_subfolder + file)

    # print(files_to_copy)

    mask_file = open(setup_subfolder + catch + "_Mask.asc", "r")

    # print(setup_subfolder + catch + "_Mask.asc")

    ncols = int(mask_file.readline().rstrip().split()[1])
    nrows = int(mask_file.readline().rstrip().split()[1])
    llx = float(mask_file.readline().rstrip().split()[1])
    lly = float(mask_file.readline().rstrip().split()[1])
    cellsize = float(mask_file.readline().rstrip().split()[1])

    mask_file.close()

    # Read in the library file from the APM simulation.
    with open(setup_subfolder + catch + "_LibraryFile.xml", 'r') as file:
        filedata = file.read()

    # Edit the library file to reflect the changed time period:
    toReplace = "<StartDay>1</StartDay>\n<StartMonth>1</StartMonth>\n<StartYear>1980</StartYear>\n<EndDay>1</EndDay" \
                ">\n<EndMonth>1</EndMonth>\n<EndYear>2011</EndYear>"

    periods = ['19801201-19901130', '19901201-20001130', '20001201-20101130', '20101201-20201130', '20201201-20301130',
               '20301201-20401130', '20401201-20501130', '20501201-20601130', '20601201-20701130', '20701201-20801130']

    replacementText = "<StartDay>01</StartDay>\n<StartMonth>12</StartMonth>\n<StartYear>1980</StartYear>\n<EndDay>30" \
                      "</EndDay>\n<EndMonth>11</EndMonth>\n<EndYear>2080</EndYear> "

    # Replace the target string
    filedata = filedata.replace(toReplace, replacementText)

    # Write the library file out again
    with open(setup_subfolder + catch + "_LibraryFile.xml", 'w') as file:
        file.write(filedata)

    # Make climate input time series files
    get_variable(
        'pr', llx, lly, cellsize, ncols, nrows, m, periods,
        setup_subfolder + "/" + catch + "_Precip.csv"
    )
    get_variable(
        'pet', llx, lly, cellsize, ncols, nrows, m, periods,
        setup_subfolder + catch + "_PET.csv"
    )
    get_variable(
        'tas', llx, lly, cellsize, ncols, nrows, m, periods,
        setup_subfolder + catch + "_Temp.csv"
    )

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
    np.savetxt(setup_subfolder + catch + "_Cells.asc", newMap, delimiter=' ', header=head, fmt='%.0f', comments='')
    # np.savetxt(newPath + "/PEAscii" + catch + ".txt", newMap, delimiter=' ', header=head, fmt='%.0f', comments='')

    return catch


# -----------------------------------------------------------------------------
# --- Create a list of catchments to process ----------------------------------
# -----------------------------------------------------------------------------

catchments = []
with open(simulation_list, 'r') as fh:
    fh.readline()
    for line in fh:
        catchments.append(line.rstrip().split(',')[0])

catchments = ["201002", "201006"]


# -----------------------------------------------------------------------------
# --- Create the necessary folder structure if needed -------------------------
# -----------------------------------------------------------------------------

if make_folder_structure:
    period_path = setup_folder + "/"
    if not os.path.exists(period_path):
        os.mkdir(period_path)
    for m in models:
        model_path = period_path + m + "/"
        if not os.path.exists(model_path):
            os.mkdir(model_path)
        for catch in catchments:
            catch_path = model_path + catch + '/'
            if not os.path.exists(catch_path):
                os.makedirs(catch_path)

# TESTING
# catch = '2002'
# time_period = 'control'
# m = '01'
# setItUp(catch, time_period, m)
# import sys
# sys.exit()

# -----------------------------------------------------------------------------
# --- Setup the catchment files -----------------------------------------------
# -----------------------------------------------------------------------------

if __name__ == '__main__':
    # print("yeah")
    # we tell python to create 4 new python instances. We can therefore run 4 tasks in parallel.
    p = Pool(processes=n_processes)  # 15
    # this python script consumes a process itself, but it's not CPU heavy so no need to factor it in.
    # A manager provides the means to safely share data with queues. Again, it consumes another process,
    # meaning we're actually going to be using 6 python instances at any one time (4 workers, this script
    # and the manager server. The manager will not consume much CPU, so again don't worry about it.
    man = Manager()
    # This is the queue, as explained earlier.
    q = man.Queue()

    # You will probably have a list of filenames to iterate over. I've just used a counter as an example

    # for time_period in ["control", "future"]:
    for m in models:
        for catch in catchments:
            # print(time_period, catch, m)

            print(catch)
            # runMe(catch, run, zipsFolder, "q")
            # we instruct workers to run runMe, using two arguments (i+1 and q). When a worker is finished,
            # workerDone will be triggered with one parameter (whatever runMe returned).

            # We're using apply_async here, because we DON'T want to WAIT before kicking off another worker.
            p.apply_async(setItUp, [catch, m, q], callback=workerDone)
            # setItUp(catch, m) # < for single processing

    # We instruct our main script to wait until all queued tasks have finished
    p.close()
    p.join()

    results = []
    while not q.empty():
        try:
            results.append(q.get())
        except:
            pass
