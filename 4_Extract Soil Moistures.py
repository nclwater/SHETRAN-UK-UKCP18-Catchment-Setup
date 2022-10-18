"""
--- Extracting Soil Moistures from SHETRAN H5 File ---

This script loads in a single H5 file using a library file and SHETRAN-IO and extracts the soil moisture data.
This data can then be fed to CEH, for use in their crop model.
You may need to install shetranio using pip - https://github.com/nclwater/shetranio

Change the historical timesteps for if you want the historical data.

Created on Thu Dec 16 11:46:40 2021
@author: Steve Birkinshaw
@editor: Ben Smith (14/10/2022)
"""


# --- PREAMBLE

import os
import numpy as np
import netCDF4 as nc
from shetranio import Model

root = "I:/SHETRAN_GB_2021/"
hist = "historical_220601_GB_APM/"
UKCP = "UKCP18rcm_220708_APM_GB/"
master_folder = root + UKCP  # hist

rcps = ['bc_01', 'bc_04', 'bc_05', 'bc_06', 'bc_07', 'bc_08',
        'bc_09', 'bc_10', 'bc_11', 'bc_12', 'bc_13', 'bc_15']

# rcps = ["bc_15"]

# Models start at 01/12/1980
spinup = 30 + 360 * 4  # 01/12/1980 + (30 + 360*4) = 01/01/1985  # Leave as 1 is you want no spinup
final_timestep = 30 + 360 * 29  # This should end in 30/12/1999


# --- DEFINE FUNCTIONS
def read_time_file(tim_file_in=str):
    """
    This will read the time file, which is useful for working out if the model is complete.
    """
    with open(tim_file_in, 'r') as time_file:
        time_line = time_file.readline().split()
        if time_line[0] == "Current":
            t = round(float(time_line[3]), 0)
        else:
            t = 0
    return t


# --- SET FOLDERS/FILES

# Run through the RCPs:
for rcp in rcps:

    rcp_no = "bc" + rcp[3:]

    # List all folders:
    folders = os.listdir(master_folder + rcp + "/")

    # Remove files from folder list:
    folders = [folder for folder in folders if "." not in folder]

    # Remove folders that have already been processed:
    nc_files = os.listdir(root + 'Simulated_Soil_Moisture/')
    nc_files = [nc.split("_")[0:2] for nc in nc_files]
    nc_files = ['_'.join(nc_f) for nc_f in nc_files]
    folders = [folder for folder in folders if rcp_no+"_"+folder not in nc_files]

    # --- PROCESS THE OUTPUTS

    # Run through each of the catchment folders:
    for folder in folders:
        files = os.listdir(master_folder + rcp + "/" + folder)

        # Find the time file:
        time_file = [file for file in files if file.endswith("tim.txt")]

        # If the simulation has not started, skip:
        if len(time_file) == 0:
            print(rcp + " - " + folder, "- not run. Skipping")
            continue

        # Else read the final timestep of the simulation:
        time = read_time_file(master_folder + rcp + "/" + folder + "/" + time_file[0])

        # If the model is incomplete, skip:
        if time < 271728:  # length of NI simulations. GB last 271752.
            print(rcp + " - " + folder, "- incomplete. Skipping")
            continue

        print(rcp + " - " + folder, "- Processing")

        # Add exception so that we can run through most of the files:
        try:

            # Find the library file:
            f = [file for file in files if 'LibraryFile.xml' in file]

            # Set the model library file to be read via SHETRAN-IO:
            model = Model(master_folder + rcp + "/" + folder + "/" + f[0])

            # # Extract the DEM for the model:
            # demvalue = model.dem

            # Open the file for writing data to (as NetCDF):
            fn = root + 'Simulated_Soil_Moisture/' + rcp_no + "_" + folder + '_sm_010185-301299.nc'
            ds = nc.Dataset(fn, 'w', format='NETCDF4')

            # Get the dimensions od the dataset:
            nrows, ncols, nlayers, _ = model.hdf.soil_moisture.values.shape
            # Set _ to final_timestep if you want to process until the end of the data.

            # Set up dataset variables:
            units = 'days since ' + str(model.hdf.soil_moisture.times[spinup - 1])
            time1 = ds.createDimension('time', final_timestep - 1 - spinup)
            ycoord = ds.createDimension('y_coord', nrows)
            xcoord = ds.createDimension('x_coord', ncols)
            times1 = ds.createVariable(units, 'f4', ('time',))
            xc = ds.createVariable('x_coord', 'f4', ('x_coord',))
            yc = ds.createVariable('y_coord', 'f4', ('y_coord',))
            value = ds.createVariable('Soil Moisture over the column', 'f4', ('y_coord', 'x_coord', 'time',))
            value.units = 'm'

            # Fill spatial variables:
            yc[:] = model.dem.y_coordinates[:]
            xc[:] = model.dem.x_coordinates[:]
            calendar = 'standard'
            for k in range(spinup, final_timestep - 1):
                date = model.hdf.soil_moisture.times[k + 1]
                times1[k - spinup] = nc.date2num(date, units, calendar)

            # set values outside catchment to zero. Otherwise the sum (-1 * -1) becomes positive
            model.hdf.vertical_thickness.square[model.hdf.vertical_thickness.square < 0] = 0

            sumOverColumn = np.zeros((nrows, ncols, final_timestep - spinup - 1))

            #  To get the soil moisture, multiply the value by the thickness (the top 5 cells down to 1.6m.
            #  0:5 is cells 0,1,2,3 and 4)
            # 'product' has dimensions (x,y,depth)
            # 'productsum' has dimensions (x,y)
            # 'sumovercolumn' has dimensions (x,y,t)

            for l in range(spinup, final_timestep - 1):

                product = np.multiply((model.hdf.vertical_thickness.square[:, :, 0:15]),
                                      (model.hdf.soil_moisture.values[:, :, 0:15, l + 1]))
                #  sum over the depths
                productsum = np.sum(product, axis=2)
                sumOverColumn[:, :, l - spinup] = productsum

            # value is total soil water going to netcdf file. Moved value out of loop to speed it up
            value[:, :, :] = sumOverColumn[:, :, :]
            ds.close()

        except Exception as e:
            print(" ---  ERROR --- ", rcp, "--" + folder, "---", e)

# ------------------------


# # --- View the Data
# import matplotlib.pyplot as plt
#
# fn = root + 'Historical_Soil_Moisture/18017_sm_010185-301299.nc'
# ds = nc.Dataset(fn, 'r', format='NETCDF4')
# ds["x_coord"][:]
# ds["y_coord"][:]
#
# print(ds)
#
# for dim in ds.dimensions.values():
#     print(dim)
#
# for var in ds.variables.values():
#     print(var)
#
# plt.imshow(ds['Soil Moisture over the column'][:, :, 10])
# plt.show()
#
# ds.close()
