# UKCP18 (12 km) Climate Change Simulations
# Ben Smith
# June/October 2022

To setup and create these simulations run:
    - 0_pet.py - to create the necessary PET data.
    - 1a_Setup... - To create non-bias-corrected simulations.
    - 2_bias_correction.... - to create new simulations of bias corrected climate data.
        **NOT USED IN FINAL OPENCLIM STUDY - see file for details**
    - 3_SHETRAN_run_simulations v2 - with copy.py - to copy files to the blades, run simulation and copy back.
    - Or, SHETRAN_copy.py & SHETRAN_run_simulations v1, which does similar.

These have the multiprocessing working - the historical probably did, but I couldn't get it working and so may have messed it up.
It works here, but you need to make sure that you're getting the right multiprocessing package - type 'python' at the start of the code
to ensure that its using the multiprocessing from your environment, not some old ArcGIS 2.7 version. I am struggling to get 1b working,
it is the multiprocessing bit that's the issue.

Further instructions are in the Setup README for UK Historical AMP, but in brief:
    - Check that you have the correct catchments listed in the Simulation_Setup_List.csv
    - Open a blade
    - Open Anaconda 3 Prompt
    - Activate conda environment: e.g. conda activate C:\ProgramData\Water_Blade_Programs\BenSmith\SHETRAN
    - Change to CONVEX folder: e.g. I:
    - Change to this directory: e.g. cd I:\SHETRAN_GB_2021\scripts\UKCP18rcm_220708_APM_GB
    - Run 1a_... : e.g. python 1a.....py - You can do this as a loop, or just set off individual RCPs separately.
    - Run 2_... : python 2_bias...py
    - Ensure that you have the run_estimates.csv how you want it so that the groupings are right.
    - Run 3_... for the simulations for your desired group: e.g. python 3_SHETRAN_run_simulations v2 - with copy.py 1
      [This copies from the convex folder to the blade - if you use V1, you need to copy manually, which may fill the blade].

    These steps may each take a long time, and you may have multiple iterations of the run script. When models crash they
    do not always get replaced. So running them, waiting a day, checking the prompt (pressing enter to clear broken models,
    then restarting the prompt and commands if necessary to reactivate up to the full number of processes. You can use
    BladeSweeper.py to sweep completed models from the blades if they have not been cleared automatically.
    
