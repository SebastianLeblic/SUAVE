# fuel_specific_excess_power.py
# 
# Created:  Jan 2023, S. Leblic (AEROCore Labs)
# 
# Modified:
# 
#           
#           
#           

"""
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

# Numpy is use extensively throughout SUAVE
import numpy as np

# Post processing plotting tools are imported here
import matplotlib.pyplot as plt
from scipy.ndimage.filters import gaussian_filter

# More basic SUAVE function
from SUAVE.Core import Data

from SUAVE.AEROCore_Additional.aero_forces_steady_state import aero_forces_steady_state
from SUAVE.AEROCore_Additional.aero_forces_steady_state import construct_aero_map
from SUAVE.AEROCore_Additional.calc_propulsion_required import calc_propulsion_required


def fuel_specific_excess_power_range(vehicle, input_details, alt_range, mach_range, aero_map_const=False):

    # unpack
    alt_vec = alt_range
    mach_vec = mach_range

    # initializing
    FPs_details = Data()
    FPs_details.FP_s = []
    FPs_details.aoa = []
    FPs_details.thrust = []
    FPs_details.TSFC = []
    FPs_details.fuel_flow_rate = []
    FPs_details.E_w = []

    # constructing aero map
    if aero_map_const:
        print("\nAERO_MAP_CONST")
        input_details.aero_map = construct_aero_map(input_details)

    # FP_s is fuel specific power
    FP_s = np.zeros(len(mach_vec))
    aoa = np.zeros(len(mach_vec))
    thrust = np.zeros(len(mach_vec))
    TSFC = np.zeros(len(mach_vec))
    fuel_flow_rate = np.zeros(len(mach_vec))
    # E_w is specific energy
    E_w = np.zeros(len(mach_vec))
    mass = 0
    speed = 0
    gravity = 0

    # iterate through range of altitude and mach to get specific power
    for i, alt in enumerate( alt_vec ):

        input_details.altitude = alt

        for j, mach in enumerate( mach_vec ):
            
            input_details.mach = mach
            
            # intake compression effects (SUAVE equations work above Mach = 0.44)
            if mach < 0.44:
                vehicle.networks.turbojet_small.inlet_nozzle.compressibility_effects = False
            elif mach >= 0.44:
                vehicle.networks.turbojet_small.inlet_nozzle.compressibility_effects = True

            # calling steady state aero/thrust force balance calculation
            aero_force_results = aero_forces_steady_state(input_details, vehicle)

            # if force-balance converged properly, aero_force_set = 0 then propulsion parameters are calculated to meet thrust required from aero-thrust force balance.
            if aero_force_results.aero_force_set == 0:
                thrust_requirements = Data()
                thrust_requirements.throttle = aero_force_results.throttle
                thrust_requirements.thrust = aero_force_results.thrust_total
                thrust_requirements.aoa = aero_force_results.aoa

                propulsion_results = calc_propulsion_required(thrust_requirements, vehicle, input_details)
                max_thrust_available = propulsion_results.max_thrust
                speed = aero_force_results.speed
                gravity = aero_force_results.gravity
                mass = vehicle.mass_properties.takeoff


                aoa[j] = aero_force_results.aoa
                thrust[j] = propulsion_results.thrust
                TSFC[j] = propulsion_results.TSFC
                fuel_flow_rate[j] = propulsion_results.fuel_flow_rate
                FP_s[j] = ((max_thrust_available - aero_force_results.drag) * speed) / ((mass * gravity) * TSFC[j])
                
                if FP_s[j] < 0:
                    FP_s[j] = 0
                    aoa[j] = 0
                    thrust[j] = 0
                    TSFC[j] = 0
                    fuel_flow_rate[j] = 0

            else:

                FP_s[j] = 0
                aoa[j] = 0
                thrust[j] = 0
                TSFC[j] = 0
                fuel_flow_rate[j] = 0
                speed = aero_force_results.speed

            # Specific energy is calculated
            PE = mass * gravity * alt
            KE = 0.5 * mass * (speed**2)
            E_w[j] = (PE + KE) / (mass * gravity)


        FPs_details.FP_s.append(FP_s.copy())
        FPs_details.aoa.append(aoa.copy())
        FPs_details.thrust.append(thrust.copy())
        FPs_details.TSFC.append(TSFC.copy())
        FPs_details.fuel_flow_rate.append(fuel_flow_rate.copy())
        FPs_details.E_w.append(E_w.copy())
    
    return FPs_details

def fuel_specific_excess_power_single_alt(vehicle, input_details, alt, mach_range):
    
    # unpack
    mach_vec = mach_range
    input_details.altitude = alt
    
    # initializing
    FPs_details = Data()
    FPs_details.FP_s = []
    FPs_details.aoa = []
    FPs_details.thrust = []
    FPs_details.TSFC = []
    FPs_details.fuel_flow_rate = []
    FPs_details.E_w = []

    # constructing aero map
    input_details.aero_map = construct_aero_map(input_details)

    # P_s is specific power
    FP_s = np.zeros(len(mach_vec))
    aoa = np.zeros(len(mach_vec))
    thrust = np.zeros(len(mach_vec))
    TSFC = np.zeros(len(mach_vec))
    fuel_flow_rate = np.zeros(len(mach_vec))
    # E_w is specific energy
    E_w = np.zeros(len(mach_vec))
    mass = 0
    speed = 0
    gravity = 0

    # iterate through range of mach to get specific power
    for j, mach in enumerate( mach_vec ):
        
        input_details.mach = mach
        
        # intake compression effects (SUAVE equations work above Mach = 0.44)
        if mach < 0.44:
            vehicle.networks.turbojet_small.inlet_nozzle.compressibility_effects = False
        elif mach >= 0.44:
            vehicle.networks.turbojet_small.inlet_nozzle.compressibility_effects = True

        # calling steady state aero/thrust force balance calculation
        aero_force_results = aero_forces_steady_state(input_details, vehicle)

        # if force-balance converged properly, aero_force_set = 0 then propulsion parameters are calculated to meet thrust required from aero-thrust force balance.
        if aero_force_results.aero_force_set == 0:
            thrust_requirements = Data()
            thrust_requirements.throttle = aero_force_results.throttle
            thrust_requirements.thrust = aero_force_results.thrust_total
            thrust_requirements.aoa = aero_force_results.aoa

            propulsion_results = calc_propulsion_required(thrust_requirements, vehicle, input_details)
            max_thrust_available = propulsion_results.max_thrust
            speed = aero_force_results.speed
            gravity = aero_force_results.gravity
            mass = vehicle.mass_properties.takeoff        
            
            aoa[j] = aero_force_results.aoa
            thrust[j] = propulsion_results.thrust
            TSFC[j] = propulsion_results.TSFC
            fuel_flow_rate[j] = propulsion_results.fuel_flow_rate
            FP_s[j] = ((max_thrust_available - aero_force_results.drag) * speed) / ((mass * gravity) * TSFC[j])

            if FP_s[j] < 0:
                FP_s[j] = 0
                aoa[j] = 0
                thrust[j] = 0
                TSFC[j] = 0
                fuel_flow_rate[j] = 0

        else:

            FP_s[j] = 0
            aoa[j] = 0
            thrust[j] = 0
            TSFC[j] = 0
            fuel_flow_rate[j] = 0
            speed = aero_force_results.speed

        # Specific energy is calculated
        PE = mass * gravity * alt
        KE = 0.5 * mass * (speed**2)
        E_w[j] = (PE + KE) / (mass * gravity)


    FPs_details.FP_s.append(FP_s.copy())
    FPs_details.aoa.append(aoa.copy())
    FPs_details.thrust.append(thrust.copy())
    FPs_details.TSFC.append(TSFC.copy())
    FPs_details.fuel_flow_rate.append(fuel_flow_rate.copy())
    FPs_details.E_w.append(E_w.copy())
    
    return FPs_details

def fuel_specific_excess_power_single_point(vehicle, input_details, alt, mach, aero_map_const=False):
    
    # unpack
    input_details.altitude = alt 
    input_details.mach = mach

    # initializing
    FPs_details = Data()
    FPs_details.FP_s = 0
    FPs_details.aoa = 0
    FPs_details.thrust = 0
    FPs_details.TSFC = 0
    FPs_details.fuel_flow_rate = 0
    FPs_details.E_w = 0
    
    # constructing aero map
    if aero_map_const:
        input_details.aero_map = construct_aero_map(input_details)

    # FP_s is fuel specific power
    FP_s = 0
    aoa = 0
    thrust = 0
    TSFC = 0
    fuel_flow_rate = 0
    # E_w is specific energy
    E_w = 0
    mass = 0
    speed = 0
    gravity = 0
    
    # intake compression effects (SUAVE equations work above Mach = 0.44)
    if mach < 0.44:
        vehicle.networks.turbojet_small.inlet_nozzle.compressibility_effects = False
    elif mach >= 0.44:
        vehicle.networks.turbojet_small.inlet_nozzle.compressibility_effects = True

    # calling steady state aero/thrust force balance calculation
    aero_force_results = aero_forces_steady_state(input_details, vehicle)

    # if force-balance converged properly, aero_force_set = 0 then propulsion parameters are calculated to meet thrust required from aero-thrust force balance.
    if aero_force_results.aero_force_set == 0:
        thrust_requirements = Data()
        thrust_requirements.throttle = aero_force_results.throttle
        thrust_requirements.thrust = aero_force_results.thrust_total
        thrust_requirements.aoa = aero_force_results.aoa

        propulsion_results = calc_propulsion_required(thrust_requirements, vehicle, input_details)
        max_thrust_available = propulsion_results.max_thrust
        speed = aero_force_results.speed
        gravity = aero_force_results.gravity
        mass = vehicle.mass_properties.takeoff
        
        aoa = aero_force_results.aoa
        thrust = propulsion_results.thrust
        TSFC = propulsion_results.TSFC
        fuel_flow_rate = propulsion_results.fuel_flow_rate
        FP_s = ((max_thrust_available - aero_force_results.drag) * speed) / ((mass * gravity) * TSFC)

        if FP_s < 0:
            FP_s = 0
            aoa = 0
            thrust = 0
            TSFC = 0
            fuel_flow_rate = 0

    else:

        FP_s = 0
        aoa = 0
        thrust = 0
        TSFC = 0
        fuel_flow_rate = 0
        speed = aero_force_results.speed

    # Specific energy is calculated
    PE = mass * gravity * alt
    KE = 0.5 * mass * (speed**2)
    E_w = (PE + KE) / (mass * gravity)


    FPs_details.FP_s = FP_s
    FPs_details.aoa = aoa
    FPs_details.thrust = thrust
    FPs_details.TSFC = TSFC
    FPs_details.fuel_flow_rate = fuel_flow_rate
    FPs_details.E_w = E_w
    
    return FPs_details