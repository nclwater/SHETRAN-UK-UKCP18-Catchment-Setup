"""
-----------------------------------------------------------------------------
Bias Correct the UKCP18 Climate Data
Ben Smith (adapted from David Prichard)
12/07/2022
-----------------------------------------------------------------------------
This script adapts the previously setup SHETRAN simulations to bias corrected
climate data. This is slightly different to the previous version on this
script as it does not break down the scenarios into decades.
Because of this, in the function for bias correcting, we no longer use the
Three distinct periods of Historical Control and Future, but instead take the
first period of the simulation that is equivalent to the historical period
and use that to replace the control. Essentially, we need a historical,
an equivalent to compare it to, and then the data we are editing. If you
change this script then make sure that you edit that second argument in the
bias correction function to ensure that it matches the historical period!

The historical period is 01/01/1980 - 01/01/2011.
The 12kmUKCP18 period is 01/12/1980 - 01/12/2080.
A difference of 335 days and a length of 10988 from the start of the UKCP18 dataset, used later.
    from datetime import date
    diff = date(1980, 12, 1) - date(1980, 1, 1)
    print(diff.days)
    print((date(2011, 1, 1) - date(1980, 12, 1)).days)

I have also changed the folder naming so that it is bc_xx, rather than xx_bc, so that the
folders are easier to sort.

This script checks what simulations have already been created via the log_bc.csv. Any that
already exist will not be reprocessed.

It is a bit weird that it runs through variables first, then catchments, but this works. Consider
changing to run through catchments, then variables.
-----------------------------------------------------------------------------
"""

import os
import shutil
import numpy as np
import pandas as pd
from scipy.stats.mstats import mquantiles
from scipy.interpolate import interp1d

# import sys

# -----------------------------------------------------------------------------

root = 'I:/SHETRAN_GB_2021/'
historical_folder = root + 'historical_220601_GB_APM/'
scenario_folder = root + 'UKCP18rcm_220708_APM_GB_NI/'
simulation_list = root + 'scripts/UKCP18rcm_220708_APM_UK/Simulation_Setup_List.csv'

log_path = scenario_folder + 'bc_log.csv'
log_append = False

# scenarios = ['control', 'future'] # [scenario]
models = ["05"]  # ['04', '06', '07', '08', '09', '10', '11', '12', '13', '15']  # '01' '05',
variables = ['Precip', 'PET', 'Temp']
variables_standard = {'Precip': 'pr', 'PET': 'pet', 'Temp': 'tas'}

use_multiprocessing = False  # True/False
nprocs = 30  # /32


# -----------------------------------------------------------------------------
def read_ascii_raster(file_path, data_type=int, return_metadata=True):
    """Read ascii raster into numpy array, optionally returning headers."""
    headers = []
    dc = {}
    with open(file_path, 'r') as fh:
        for i in range(6):
            line = fh.readline()
            headers.append(line.rstrip())
            key, val = line.rstrip().split()
            dc[key.lower()] = val
    ncols = int(dc['ncols'])
    nrows = int(dc['nrows'])
    xll = float(dc['xllcorner'])
    yll = float(dc['yllcorner'])
    cellsize = float(dc['cellsize'])
    nodata = float(dc['nodata_value'])  # Changed from NODATA_value

    arr = np.loadtxt(file_path, dtype=data_type, skiprows=6)

    headers = '\n'.join(headers)
    headers = headers.rstrip()

    if return_metadata:
        return arr, ncols, nrows, xll, yll, cellsize, headers
    else:
        return arr


def eqm(obs, p, s, nbins=10, extrapolate=None):
    """Empirical quantile mapping.

    Based on: https://svn.oss.deltares.nl/repos/openearthtools/trunk/python/applications/hydrotools/hydrotools/statistics/bias_correction.py

    Args:
        obs: observed climate data for the training period
        p: simulated climate by the model for the same variable obs for the
            training ("control") period
        s: simulated climate for the variables used in p, but considering the
            test/projection ("scenario") period
        nbins: number of quantile bins
        extrapolate: None or 'constant', indicating the extrapolation method to
            be applied to correct values in 's' that are out of the range of
            lowest and highest quantile of 'p'

    """
    binmid = np.arange((1. / nbins) * 0.5, 1., 1. / nbins)

    qo = mquantiles(obs[np.isfinite(obs)], prob=binmid)
    qp = mquantiles(p[np.isfinite(p)], prob=binmid)

    p2o = interp1d(qp, qo, kind='linear', bounds_error=False)

    c = p2o(s)

    if extrapolate is None:
        c[s > np.max(qp)] = qo[-1]
        c[s < np.min(qp)] = qo[0]
    elif extrapolate == 'constant':
        c[s > np.max(qp)] = s[s > np.max(qp)] + qo[-1] - qp[-1]
        c[s < np.min(qp)] = s[s < np.min(qp)] + qo[0] - qp[0]

    return c


def process_catchment(catch, model, variable):  # , q=None
    print(model, variables_standard[variable], catch)
    # --- Copy pre-correction files that do not need to be changed
    destination_folder = scenario_folder + "bc_" + model + '/' + catch + "/"
    # E.g. "I:/SHETRAN_GB_2021/UKCP18rcm_220708_APM_GB/01_bc/1001/"

    if not os.path.exists(destination_folder):
        os.mkdir(destination_folder)
        # TODO Make this work recursively, currently it doesn't work unless the bc_xx is already created.

    if variable == 'Precip':

        # - cells map needs to come from historical, not pre-correction
        src = historical_folder + catch + '/' + catch + '_Cells.asc'
        dst = destination_folder + catch + '_Cells.asc'
        shutil.copy(src, dst)

        # - rest of files
        source_folder = scenario_folder + model + '/' + catch + '/'
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
            shutil.copy(source_folder + file, destination_folder + file)

    # Read catchment mask
    mask_path = historical_folder + catch + '/' + catch + '_Mask.asc'
    mask, ncols, nrows, xll, yll, cellsize, hdrs = read_ascii_raster(
        mask_path, data_type=int, return_metadata=True
    )

    # Read historical and scenario cell ID files
    historical_cells_path = historical_folder + catch + '/' + catch + '_Cells.asc'
    scenario_cells_path = scenario_folder + model + '/' + catch + '/' + catch + '_Cells.asc'

    historical_cells = read_ascii_raster(historical_cells_path, return_metadata=False)
    scenario_cells = read_ascii_raster(scenario_cells_path, return_metadata=False)

    # Read UKCP18 pre-correction time series file:
    subfolder = scenario_folder + model + '/' + catch + '/'
    uncorr_ukcp18_series = pd.read_csv(subfolder + catch + '_' + variable + '.csv')

    # Read historical time series file
    historical_subfolder = historical_folder + catch + '/'
    historical_path = historical_subfolder + catch + '_' + variable + '.csv'
    historical_series = pd.read_csv(historical_path)

    # Cell-wise quantile mapping
    corr_series = pd.DataFrame()
    for yi in range(nrows):
        for xi in range(ncols):
            if mask[yi, xi] == 0:
                h_cell = str(historical_cells[yi, xi])
                s_cell = str(scenario_cells[yi, xi])

                h_series = historical_series.loc[:, h_cell].values
                uc_series = uncorr_ukcp18_series.loc[:, s_cell].values

                c_series = eqm(obs=h_series,  # Observed data
                               p=uc_series[335:10988],  # Training period (01/11/1980-31/12/2010?)
                               s=uc_series,  # Scenario period (i.e. whole UKCP18 period)
                               nbins=100,  # I decreased bins. I don't think it makes any real difference.
                               extrapolate='constant')

                # Make sure that we have not bias controlled any negative values:
                if variable in ['Precip', 'PET']:
                    c_series[c_series < 0.0] = 0.0

                corr_series = pd.concat((corr_series, pd.Series(c_series, name=h_cell)), axis=1)
                # TODO Check that this works if the columns are not in the correct order...

    # Write the corrected dataset:
    corr_series = corr_series.round(2)
    output_path = destination_folder + catch + '_' + variable + '.csv'
    corr_series.to_csv(output_path, index=False)


def process_catchment_mp(catch, model, variable, q):
    try:
        process_catchment(catch, model, variable)
        output_line = [model, variable, catch, 'Y']
    except Exception as e:
        print(e)
        output_line = [model, variable, catch, 'N']
    output_line = ','.join(output_line)
    q.put(output_line)


def log_status(log_path, q):
    with open(log_path, 'a') as fh:
        if not log_append:
            fh.write('Model,Variable,Catchment,Flag\n')
        while True:
            msg = q.get()
            if msg == 'kill':
                break
            fh.write(msg + '\n')
            fh.flush()


def process_mp(cases_to_processs):
    manager = mp.Manager()
    q = manager.Queue()
    pool = mp.Pool(nprocs)

    logger = pool.apply_async(log_status, (log_path, q))

    jobs = []

    for case in cases_to_process:
        model, variable, catch = case
        job = pool.apply_async(process_catchment_mp, (catch, model, variable, q))
        jobs.append(job)

    for job in jobs:
        job.get()

    q.put('kill')
    pool.close()
    pool.join()


# -----------------------------------------------------------------------------

catchments = []
with open(simulation_list, 'r') as fh:
    fh.readline()
    for line in fh:
        catchments.append(line.rstrip().split(',')[0])
# catchments=[""]

# Make initial log file based on cases processed so far
# - THEN COMMENT OUT
## models = ['01', '04', '05']
## with open(log_path, 'a') as fh_log:
##     fh_log.write('Scenario,Model,Variable,Catchment,Flag\n')
##     for scenario in scenarios:
##         for model in models:
##             for variable in variables:
##                 for catch in catchments:
##                     print(scenario, model, variable, catch)
##                     case_folder = scenario_folder + '/' + scenario + '_bc/' + model + '/' + catch + '/'
##                     series_path = case_folder + catch + '_' + variable + '.csv'
##                     if os.path.exists(series_path):
##                         flag = 'Y'
##                     else:
##                         flag = 'N'
##                     output_line = [scenario, model, variable, catch, flag]
##                     output_line = ','.join(output_line)
##                     fh_log.write(output_line + '\n')
## sys.exit()

processed = []
if log_append:
    with open(log_path, 'r') as fh_log:
        fh_log.readline()
        for line in fh_log:
            line = line.rstrip().split(',')
            ##case = line[:-1]
            case = line
            case = '_'.join(case)
            processed.append(case)

cases_to_process = []
for model in models:
    for catch in catchments:
        for variable in variables:
            case = '_'.join([model, variable, catch, 'Y'])
            if case not in processed:
                cases_to_process.append([model, variable, catch])

if __name__ == "__main__":

    if use_multiprocessing:
        import multiprocessing as mp

        process_mp(cases_to_process)
    else:
        with open(log_path, 'a') as fh_log:
            if not log_append:
                fh_log.write('Model,Variable,Catchment,Flag\n')

        for case in cases_to_process:
            model, variable, catch = case

            try:
                process_catchment(catch, model, variable)
                output_line = [model, variable, catch, 'Y']
                output_line = ','.join(output_line)
                fh_log.write(output_line + '\n')
            except Exception as e:
                print(e)
                output_line = [model, variable, catch, 'N']
                output_line = ','.join(output_line)
                fh_log.write(output_line + '\n')
