# max_mach_at_alt.py
#
# Created: Jan 2023, S. Leblic
#
# Modified:
#
#
#
#

"""

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

from SUAVE.AEROCore_Additional              import aero_forces
from SUAVE.AEROCore_Additional.aero_forces  import construct_aero_map
from SUAVE.AEROCore_Additional              import calc_propulsion_required
from SUAVE.Core                             import Data

import numpy as np


def max_mach_at_alt(input_details, vehicle):
    
    # settings
    error = 1e-7

    #initialize 
    mach_high = 1.3
    mach_low = 0.8
    mach_delta = 50

    # constructing aero map
    input_details.aero_map = construct_aero_map(input_details)

    while mach_delta > error:

        mach_middle = (mach_high + mach_low) / 2

        # intake compression effects (SUAVE equations work above Mach = 0.44)
        if mach_middle < 0.44:
            vehicle.networks.turbojet_small.inlet_nozzle.compressibility_effects = False
        elif mach_middle >= 0.44:
            vehicle.networks.turbojet_small.inlet_nozzle.compressibility_effects = True

        input_details.mach = mach_middle
        aero_results = aero_forces(input_details, vehicle)
        
        thrust_requirements = Data()
        thrust_requirements.throttle = aero_results.throttle
        thrust_requirements.thrust = aero_results.thrust_total
        thrust_requirements.aoa = aero_results.aoa

        propulsion_results = calc_propulsion_required(thrust_requirements, vehicle, input_details)
        throttle = propulsion_results.throttle

        if throttle > 1:
            mach_high = mach_middle
        elif throttle < 1:
            mach_low = mach_middle
        
        mach_delta = np.absolute(mach_high - mach_low)

    max_mach = mach_middle
    
    return max_mach, propulsion_results, aero_results
        