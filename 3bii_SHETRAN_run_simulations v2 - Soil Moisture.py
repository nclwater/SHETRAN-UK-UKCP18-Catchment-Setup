"""
Batch Run SHETRAN Simulations on the Blades, but storing them elsewhere
Ben Smith (adapted from David Prichard)
18/10/2022

This script is based on 3_SHETRAN_run_simulations v2. It sets up the catchments to run for only the baseline period and
outputs are reduced to only soil moisture. This is because the previous setup originally only output for the top 6
layers, which wasn't enough (it was designed using the autocalibration version, which has thicker layers).

Notes:
    - this is taken from:
        CONVEX/SHETRAN_GB_2021/scripts/historical_220601_SHETRAN_UK/Run Scripts
    - Catchments are split into groups in the runtimes_path file. This is so that
        a different group can run on each blade. This script is run from the command line
        and the first argument (a number) is taken as the group to process when running this script.
        E.g. >> python runs.py 1 # sets group_to_process to 1
    - I normally run this script from the CONVEX Scripts folder.
    - This has been edited, so now it will, hopefully, run through the different scenarios in
        turn, as specified in the simulation_folder LIST.
    - This will also edit the visualisation plan in a different way to before. The chaps
        from CEH wanted daily soil moisture for the UK for the baseline period (about 1980-2010),
        they also only want these from the UKCP18 runs, so we only need a limited period of recording.
        As such, we will add another TIME to the visualisation plan and limit it to 30 years of output.
    - The runtime estimates here are split into two groups only. This is because they will multiprocess on the blades.

UPDATES in V2
    - The blades are quite full, and so runs crash when disk space runs out. Therefore, this script copies files from
        CONVEX as they are needed, runs the simulation, then copies back to convex and deletes the files on the blades.
            (Mainly an issue on Blade 1)
    - If there are already files in the source simulation folder, these WILL be overwritten. This means that
        simulations WILL be updated if rerun. You can change this by editing the 'overwrite' variable.

"""
import os
import shutil
import subprocess
import sys
import time
import pandas as pd

# -----------------------------------------------------------------------------
# USER INPUTS

# Set file paths:
executable_folder = "C:/ProgramData/Water_Blade_Programs/BenSmith/SHETRAN_snow"

CONVEX_simulation_folder = "I:/SHETRAN_GB_2021/05_Climate_Change_Simulations/UKCP18rcm_181022_APM_UK_CEH/"
blade_folder = "C:/BenSmith/Blade_SHETRANGB_OpenCLIM_UKCP18rcm_220708_APM/Temp_simulations/CEH/"

if not os.path.isdir(blade_folder):
    os.mkdir(blade_folder)

# Set climate runs to simulate:
# rcm_suites = ["bc_04", "bc_05", "bc_06", "bc_07", "bc_08", "bc_09", "bc_10", "bc_11", "bc_12", "bc_13", "bc_15"]
rcm_suites = ["bc_13"]

# Set the file containing catchment names and multiprocessing groups:
runtimes_path = "I:/SHETRAN_GB_2021/05_Climate_Change_Simulations/UKCP18rcm_181022_APM_UK_CEH/Catchment Ranking for " \
                "Soil Moisture.csv"

# Setup input argument for taking the group that you would like to process from conda prompt:
group_to_process = int(sys.argv[1])
# ^^ Type in the number of the group in the runtime estimates that you wish to process in the
# conda prompt after the name of this script.

# Choose whether you want multiprocessing:
use_multiprocessing = True
num_processes = 12

# Read in the catchment/simulation name and their processing groups:
catchments = pd.read_csv(runtimes_path)
catchments = catchments[catchments["Run_Group"] == group_to_process]
catchments = catchments["Catchment"].to_list()
catchments = [str(cat) for cat in catchments]


# For testing:
# catchments = ["54022"]

# -----------------------------------------------------------------------------
# FUNCTIONS FOR RUNNING SIMULATIONS

def run_SHETRAN_batch_with_copy(catchment, simulation_run_folder, simulation_store_folder=None, q=None):
    print("----------------------------------------------------------------------------------")
    print("SIMULATION: " + catchment)
    print("----------------------------------------------------------------------------------")

    # Copy simulation to Blade (creating directory) if doing that:
    if simulation_store_folder is not None:
        folder_copy(simulation_store_folder, simulation_run_folder, True)

    os.chdir(str(executable_folder))
    start_time = time.time()

    # Run the prepare script:
    try:
        return_code = subprocess.call(
            ['./shetran-prepare-snow.exe', simulation_run_folder + catchment + "_LibraryFile.xml"])
        prep_flag = 'Y'
    except Exception as e:
        print(e)
        prep_flag = 'N'

    # Edit the Visualisation plan to get daily soil moisture in the top 2m for 30 years:
    visualisation_plan_swap_line(
        "GRID_OR_LIST_NO^7 : TIMES^8 : LAYERS^1 1 : ENDITEM",
        "GRID_OR_LIST_NO^7 : TIMES^10 : LAYERS^1 15 : ENDITEM",
        simulation_run_folder + "input_" + catchment + "_visualisation_plan.txt")

    visualisation_plan_swap_line(
        "stop",
        "times\n10 1 !number and no. of entries\n24 260064 !every 24 hours for 30.1 years (1980-2010)\n\nstop",
        simulation_run_folder + "input_" + catchment + "_visualisation_plan.txt")

    # Edit the Visualisation plan to remove the snow output item:
    visualisation_plan_remove_item("6", simulation_run_folder + "input_" + catchment + "_visualisation_plan.txt")
    visualisation_plan_remove_item("5", simulation_run_folder + "input_" + catchment + "_visualisation_plan.txt")
    visualisation_plan_remove_item("1", simulation_run_folder + "input_" + catchment + "_visualisation_plan.txt")
    # Edits to crop outputs further for CEH:
    visualisation_plan_remove_item("3", simulation_run_folder + "input_" + catchment + "_visualisation_plan.txt")
    visualisation_plan_remove_item("1", simulation_run_folder + "input_" + catchment + "_visualisation_plan.txt")

    # Run the SHETRAN simulation:
    try:
        return_code = subprocess.call(  # Requires the 3 separate arguments
            ['shetran.exe', '-f ', simulation_run_folder + "rundata_" + catchment + ".txt"])
        run_flag = 'Y'
    except:
        run_flag = 'N'

    # Log the simulation completion status:
    simulation_final_timestep = read_time_file(simulation_run_folder + "output_" + catchment + "_tim.txt")

    # Copy and delete from Blade:
    if simulation_store_folder is not None:
        folder_copy(simulation_run_folder, simulation_store_folder, overwrite=True)
        shutil.rmtree(simulation_run_folder)

    print("- - - - - - - - - - - - - - - - - - - -")
    print("COMPLETED (OR CRASHED): " + catchment)
    print("- - - - - - - - - - - - - - - - - - - -")

    # Log simulation:
    simulation_log_entry = catchment + ',' + prep_flag + ',' + run_flag + ',' + str(
        round((time.time() - start_time) / 60, 2)) + ',' + str(simulation_final_timestep)
    if use_multiprocessing:
        q.put(simulation_log_entry)
    else:
        return simulation_log_entry


def run_SHETRAN_batch_with_copy_mp(catchment_list, simulation_source_folder,
                                   temp_simulation_folder, simulation_subfolder=""):
    manager = mp.Manager()

    q = manager.Queue()

    pool = mp.Pool(num_processes)

    logger = pool.apply_async(log_progress, (q,))

    # Set subfolder as appropriate, according to whether there is a subfolder of folders of catchments / simulations.
    # if simulation_subfolder is None:
    #     simulation_subfolder = ""

    jobs = []
    for catchment in catchment_list:
        simulation_source_subfolder = simulation_source_folder + simulation_subfolder + "/" + catchment + "/"
        temp_folder = temp_simulation_folder + simulation_subfolder + "_" + catchment + "/"

        # Set up a 'job', i.e. a simulation:
        job = pool.apply_async(run_SHETRAN_batch_with_copy, (catchment, temp_folder, simulation_source_subfolder, q))
        # Arguments used above are: "catchment", "simulation_run_folder", "simulation_store_folder", "q".

        # Add that job to the list of jobs:
        jobs.append(job)

    # Run the simulations (i.e. jobs)
    try:
        for job in jobs:
            try:
                # Adding this try here means that if there is a general (i.e. non-SHETRAN) error,
                # it will continue to try getting jobs. Else it will just crash. E.g. missing library file.
                job.get()
            except:
                print("ERROR in one of the job files. Skipping. Perhaps a missing simulation file?")

    # But include the option to kill them - 1st attempt, may need refining... but seems to work
    except KeyboardInterrupt:
        print(">  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  <")
        print("> Caught KeyboardInterrupt, terminating workers <")
        print(">  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  <")
        pool.terminate()

    # q.put('kill')  # I don't know what this does... but it seems to work against the terminate above...
    pool.close()
    pool.join()


def folder_copy(source_folder, destination_folder, overwrite=False):
    if not os.path.isdir(destination_folder):
        os.mkdir(destination_folder)

    files_2_copy = os.listdir(source_folder)

    if not overwrite:
        destination_files = os.listdir(destination_folder)
        files_2_copy = [i for i in files_2_copy if i not in destination_files]

    for file in files_2_copy:
        shutil.copy2(source_folder + file, destination_folder + file)

    return files_2_copy


def visualisation_plan_swap_line(old_line, new_line, file_in, file_out=None, strip_ws=True):
    """
    old_line & new_line  - Strings of the full lines in the visualisation plan (without white space,
                            see strip_ws).
    file_out             - Do not specify  if you want to overwrite.
    strip_ws             - True/False depending on whether you want trailing white space to be matched and included
                            in output. - Default True, so do not include white space in old line (and probably not new
                             line, just for consistency).
    ALSO, consider whether there are multiple matches.

    TODO Make a replacement method based on line number instead of string matching.
    """

    if file_out is None:
        file_out = file_in

    with open(file_in, 'r') as vis:

        replacement = ""
        change_checker = 0

        for line in vis:

            if strip_ws:
                line = line.rstrip()

            changes = line.replace(old_line, new_line)

            if line != changes:
                change_checker += 1
            replacement = replacement + changes + "\n"

    with open(file_out, "w") as new_vis:
        new_vis.write(replacement)

    if change_checker == 0:
        return "WARNING: No changes made"


def visualisation_plan_remove_item(item_number, vis_file_in=str, vis_file_out=None):
    """
    Don't forget that if you use this is combination with the number altering that you need to match the altered number.
    If you are removing multiple items, remove the higher numbers first.
    item_number can be a string or integer.
    Do not specify file_out if you want to overwrite.
    """

    if vis_file_out == None:
        vis_file_out = vis_file_in

    with open(vis_file_in, 'r') as vis:
        updated_text = ""
        number_corrector = 0

        for line in vis:
            line = line.strip().split(" : ")

            # IF the line starts with item then skip ('item' will be written later)
            if line[0] == "item":
                continue

            # IF the line starts with NUMBER, decide whether to read or write:
            if line[0][0:len(line[0]) - 2] == "NUMBER":

                # IF it is the number of interest read the next line too, not writing either
                # and add one to the index corrector:
                if line[0][-1] == str(item_number):
                    next(vis)
                    number_corrector += 1

                # IF a different number:
                if line[0][-1] != str(item_number):
                    new_number = int(line[0][-1]) - number_corrector
                    line[0] = str(line[0][0:len(line[0]) - 1] + str(new_number))
                    updated_text = updated_text + 'item \n' + " : ".join(line) + "\n" + next(vis)

            # If neither, just copy the line:
            else:
                updated_text = updated_text + " : ".join(line) + "\n"

    with open(vis_file_out, "w") as new_vis:
        new_vis.write(updated_text)

    if new_number == 0:
        return "WARNING: No lines were edited"


def read_time_file(tim_file_in=str):
    with open(tim_file_in, 'r') as time_file:

        time_line = time_file.readline().split()

        if time_line[0] == "Current":
            t = str(round(float(time_line[3]), 0))
        else:
            t = 0
    return t


def log_progress(q):
    with open(log_path, 'w') as fh:
        fh.write('Catchment, prep_ReturnCode, run_ReturnCode, Simulation Duration (mins), Model Timestep (hrs) \n')
        while True:
            msg = q.get()
            if msg == 'kill':
                break
            fh.write(msg + '\n')
            fh.flush()


# -----------------------------------------------------------------------------
# CODE FOR RUNNING THROUGH THE SIMULATIONS EITHER IN PARALLEL OR IN SERIES

# Launch runs
if __name__ == '__main__':
    if use_multiprocessing:

        import multiprocessing as mp

        for rcm in rcm_suites:
            run_SHETRAN_batch_with_copy_mp(catchment_list=catchments,
                                           simulation_source_folder=CONVEX_simulation_folder,
                                           temp_simulation_folder=blade_folder,
                                           simulation_subfolder=rcm)

    else:
        # TODO I don't know whether this works with the copy as it doesn't define the temp_simulation_folder...
        for rcm in rcm_suites:
            log_path = CONVEX_simulation_folder + rcm + '/runlog_' + str(group_to_process) + '.txt'
            write_log_header = True
            for c in range(len(catchments)):
                log_entry = run_SHETRAN_batch_with_copy(
                    catchment=catchments[c],
                    simulation_run_folder=CONVEX_simulation_folder + rcm + "/" + catchments[c] + '/')
                with open(log_path, 'w') as fh:
                    if write_log_header:
                        fh.write('Catchment, prep_ReturnCode, run_ReturnCode, Simulation Duration (mins), '
                                 'Model Timestep (hrs) \n')
                        write_log_header = False
                    fh.write(log_entry + '\n')
                    fh.flush()
