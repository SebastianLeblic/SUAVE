# max_mach_at_alt_setup_template.py
#
# Created: Jan 2023, S. Leblic (AEROCore Labs)
#
# Modified:
#
#
#
#

"""

TEMPLATE FILE

Finds maximum mach for a given aircraft at a given altitude.


Inputs:
    - SUAVE vehicle
    - Aerodynamics coefficient file (or CFD file)

Outputs:
    - Max mach
    - Aerodynamics:
        - Angle of attack [degrees]
        - Drag force [N]
        - Lift force [N]

    - Propulsion:
        - Thrust [N]
        - Throttle [%]
        - Fuel flow rate [kg/s]
        - TSFC [g/kN.s]
    

"""

# ---------------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------------

import SUAVE

from SUAVE.AEROCore_Additional              import max_mach_at_alt
from SUAVE.Core                             import Data

import numpy as np

# setup vehicle
import sys, os
from os import path

# changing directory path to extract vehicle files. Uncomment print basepath statement to see current path and adjust append as necessary
basepath = path.dirname(__file__)
#print("\n basepath:", basepath)
filepath = path.abspath(path.join(basepath, "..", "Vehicles"))
sys.path.append(filepath)
from VEHICLE_NAME_HERE import vehicle_setup



def main():
    
    # settings
    alt = SET_ALTITUDE_HERE # units [m], set altitude to be assessed here.

    # initialize 
    input_details = Data()
    input_details.aero_file = basepath + '\CFD_Results\\AERO_FILE_NAME_HERE.txt'
    input_details.aero_plotting = False
    input_details.d_isa = 0
    input_details.altitude = alt

    # construct the baseline aircraft
    vehicle = full_setup()

    max_mach, propulsion_results, aero_results = max_mach_at_alt(input_details, vehicle)
    
    print("\n propulsion_results:", propulsion_results)
    print("\n aero_results:", aero_results)
    print("\n altitude (m):", alt)
    print("\n max_mach:", max_mach)
    
    return

def full_setup():
    
    # Vehicle data
    vehicle  = vehicle_setup()
       
    
    return vehicle
        
if __name__ == '__main__': 
    main()
        