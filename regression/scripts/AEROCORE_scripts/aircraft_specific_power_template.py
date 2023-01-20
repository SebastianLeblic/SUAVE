# aircaft_specific_power_template.py
# 
# Created:  Jan 2023, S. Leblic (AEROCore Labs)
# 
# Modified:           
#           
#           
#           

"""

TEMPLATE FILE

Assesses the specific excess power of a given aircraft for a range of speeds and altitudes. Inherently a performance envelope is created, where max Mach can be determined. This further acts as a first step to developing an energy based approach to optimal trajectory mapping.

Inputs:
    - SUAVE vehicle (for propulsion system description)
    - Aerodynamics coefficient file

Outputs:
    - Specific excess power map for a range of altitudes and Mach
    - Max acheivable Mach number, and at what altitude

"""

# ----------------------------------------------------------------------
#   Imports
# ----------------------------------------------------------------------

import SUAVE
# Units allow any units to be specificied with SUAVE then automatically converting them the standard
from SUAVE.Core import Units
from SUAVE.Plots.Performance.Mission_Plots import * 

# Numpy is use extensively throughout SUAVE
import numpy as np

# More basic SUAVE function
from SUAVE.Core import Data

import time

import sys, os
from os import path

# changing directory path to extract vehicle files. Uncomment print basepath statement to see current path and adjust append statement as necessary
basepath = path.dirname(__file__)
#print("\n basepath:", basepath)
filepath = path.abspath(path.join(basepath, "..", "Vehicles"))
sys.path.append(filepath)
from VEHICLE_NAME_HERE import vehicle_setup   # indicate vehicle name as SUAVE vehicle file name from "Vehicles" folder

from SUAVE.AEROCore_Additional.specific_excess_power import specific_excess_power

# ----------------------------------------------------------------------
#   Main
# ----------------------------------------------------------------------

def main():
    
    st = time.time()
    

    # construct the baseline aircraft
    vehicle = full_setup()

    # search space settings
    n = 50
    alt_range = np.linspace(0, 16000, n)
    mach_range = np.linspace(0.01, 1.2, n)
    
    # initializing and settings
    input_details = Data()
    # change CFD file name where indicated here
    input_details.aero_file = basepath + '\CFD_Results\\CFD_FILE_NAME_HERE.txt'
    input_details.aero_plotting = False     # set True if viewing aero map plots is desired
    input_details.Ps_plotting = True        # set True if specific excess power and other additional plots are desired
    input_details.Ps_plot_smoothing = False # set True if n is sparse and smoothing of contour plots is desired
    input_details.save_ps_plots = False     # set True if saving Ps_plots is desired. Saved in .pdf format in same folder location
    input_details.d_isa = 0                 # for propulsion calcs - set temperature deviation from standard temp setting (default 298.15 K)

    # call function to formulate specific excess power map
    P_s_details = specific_excess_power(vehicle, input_details, alt_range, mach_range)

    # print details of point of maximum acheivable Mach
    print("\n max_mach:", P_s_details.max_mach)
    print("\n max_mach_alt:", P_s_details.max_mach_alt)

            
    et = time.time()
    elapsed_time = et - st
    print("\n Execution time: ", elapsed_time, 's')

    
    return


# ----------------------------------------------------------------------
#   Analysis Setup
# ----------------------------------------------------------------------

def full_setup():
    
    # acquire vehicle data
    vehicle  = vehicle_setup()
       
    
    return vehicle
        
        
if __name__ == '__main__': 
    main()
        
