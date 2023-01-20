# energy_climb_trajectory_const_mass.py
# 
# Created: Jan 2023, S.Leblic (AEROCore Labs)
#
# Modified:
#
#
#
#

"""

Calculates either minimum time-to-climb or minimum fuel-to-climb using energy methods as outlined by Rutowski (1954), "Energy Approach to the General Aircraft Performance Problem".

The produced trajectory profile can be used in a mission as a path for the aircraft to follow. Aircraft maintains a constant mass throughout.

Inputs:
    - SUAVE vehicle description
    - aero coefficients map
    - starting altitude     [m]
    - starting mach         []

Outputs:
    - Climb trajectory (altitude and mach points)


"""

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------

import numpy as np
from scipy.optimize     import curve_fit

from SUAVE.AEROCore_Additional.specific_excess_power import *
from SUAVE.AEROCore_Additional.fuel_specific_excess_power import *
from SUAVE.AEROCore_Additional.specific_energy import specific_energy



def minimum_fuel_to_climb_const_mass(vehicle, input_details):
    
    # find point of maximum mach
    n = 100     # this should be set around >100 to sufficiently refine search space
    alt_range_init = np.linspace(6000, 15000, n)
    mach_range_init = np.linspace(0.8, 1.1, n)

    FPs_init_details = specific_excess_power_range(vehicle, input_details, alt_range_init, mach_range_init)
    init_alt = input_details.init_alt
    init_mach = input_details.init_mach
    mach_final = FPs_init_details.max_mach
    alt_final = FPs_init_details.max_mach_alt
    
    # find energy height of final mach and altitude
    E_w_final = specific_energy(vehicle, alt_final, mach_final)

    # at constant altitude from trajectory starting point, find point of maximum specific excess power; this is the second point of the trajectory
    n1 = 40
    mach_range_temp = np.linspace(init_mach, mach_final, n1)
    FPs_max_init = fuel_specific_excess_power_single_alt(vehicle, input_details, init_alt, mach_range_temp)
    i, j = np.where(FPs_max_init.FP_s == np.max(FPs_max_init.FP_s))
    E_w_init = specific_energy(vehicle, init_alt, init_mach)

    # initialize trajectory data structure. All trajectory points (alt, mach) are stored in this class.
    trajectory = Data()
    trajectory.mach = []
    trajectory.altitude = []
    trajectory.E_w = []

    trajectory.mach.append(init_mach)
    trajectory.mach.append(np.float(mach_range_temp[j]))
    trajectory.altitude.append(init_alt/1000)
    trajectory.altitude.append(init_alt/1000)
    trajectory.E_w.append(np.float(E_w_init))
    
    # define planet and atmospheric models for gravity and air property calculations
    planet = SUAVE.Attributes.Planets.Earth()
    atmo = SUAVE.Analyses.Atmospheric.US_Standard_1976()

    # increments
    alt_increment = 500     # units [m]
    E_w_increment = 250     # units [m] (specific energy can be called energy height)

    # energy height of second trajectory point (begin climb trajectory from here)
    E_w_current = specific_energy(vehicle, init_alt, trajectory.mach[1])
    trajectory.E_w.append(np.float(E_w_current))
    current_mach = trajectory.mach[-1]
    alt = init_alt

    # initialize toggle. Toggle is used for identifying certain portions of climb trajectory.
    E_w_set = 0

    while E_w_set != 2:

        # if aircraft specific energy is equal to that at final altitude and mach
        if E_w_set == 1:
            E_w_set = 2

        # if specific energy increment goes past that at final altitude and mach, set next energy height to final energy height
        if (E_w_final - E_w_current) < E_w_increment and E_w_set != 2:
            E_w_new = E_w_final
            E_w_set = 1

        # increase energy height by increment
        elif E_w_set == 0:
            E_w_new = E_w_current + E_w_increment
            

        if E_w_set != 2:
            # calculate contour line at new energy height
            n2 = 20
            alt_array = np.linspace((alt - (alt_increment / 2)), (alt + alt_increment + 500), n2)
            mach_array = np.zeros(n2)

            # for given energy height, calculate subsequent speed at a range of altitudes along the contour
            for r, height in enumerate(alt_array):
                gravity_new = planet.compute_gravity(height)
                V = np.sqrt((E_w_new - height) * 2 * gravity_new)
                air_properties_new = atmo.compute_values(height, 0.)
                mach_array[r] = V / air_properties_new.speed_of_sound
            
            # cubic polynomial approximation
            def objective(x, a, b, c, d):
                return a + (b * x) + (c * x**2) + (d * x**3)

            # using cubic polynomial approximation of energy height contour
            popt1, _ = curve_fit(objective, mach_array, alt_array)
            a, b, c, d = popt1
            
            # initialize fuel specific excess power array
            FP_s_array = np.zeros(len(mach_array))

            # calculate specific excess power along constant energy height contour
            for o, mach in enumerate(mach_array):
                alt_temp = alt_array[o]
                FP_s_array[o] = fuel_specific_excess_power_single_point(vehicle, input_details, alt_temp, mach).FP_s

            # cubic polynomial approximation of specific excess power curve
            popt2, _ = curve_fit(objective, mach_array, FP_s_array)
            e, f, g, h = popt2

            mach_high = mach_array[0]
            mach_low = mach_array[-1]
            mach_delta = 10

            # find point of zero slope along specific excess power curve. Using the bisection method, the mach at the point of zero slope is determined
            while mach_delta > 1e-8:
                FP_s_der_high = f + ((2 * g) * mach_high) + ((h * 3) * mach_high**2)

                mach_mid = (mach_high + mach_low) / 2

                FP_s_der_mid = f + ((2 * g) * mach_mid) + ((h * 3) * mach_mid**2)

                if FP_s_der_mid * FP_s_der_high > 0:
                    mach_high = mach_mid
                elif FP_s_der_mid * FP_s_der_high < 0:
                    mach_low = mach_mid
                
                mach_delta = mach_high - mach_low

            
            current_mach = mach_mid
            
            # calculate altitude along energy height contour for given mach
            alt = objective(current_mach, a, b, c, d)
            E_w_current = specific_energy(vehicle, alt, current_mach)
            
            trajectory.mach.append(current_mach)
            trajectory.altitude.append(alt / 1000)
            trajectory.E_w.append(np.float(E_w_current))
        
        # if final energy height is achieved, follow energy height contour to final mach and altitude
        else:
            n2 = 20
            alt_array = np.linspace(alt, alt_final, n2)
            mach_array = np.zeros(n2)

            for p, height in enumerate(alt_array):
                    gravity_new = planet.compute_gravity(height)
                    V = np.sqrt((E_w_new - height) * 2 * gravity_new)
                    air_properties_new = atmo.compute_values(height, 0.)
                    mach_array[p] = V / air_properties_new.speed_of_sound
                    trajectory.mach.append(mach_array[p])
                    trajectory.altitude.append(height / 1000)
                    trajectory.E_w.append(np.float(E_w_new))
                    E_w_current = E_w_final

        
    trajectory.mach_final = mach_final
    trajectory.alt_final = alt_final

    return trajectory

def minimum_time_to_climb_const_mass(vehicle, input_details):

    # find point of maximum mach
    n = 100     # this should be set around >100 to sufficiently refine search space
    alt_range_init = np.linspace(6000, 15000, n)
    mach_range_init = np.linspace(0.8, 1.1, n)

    Ps_init_details = specific_excess_power_range(vehicle, input_details, alt_range_init, mach_range_init)
    init_alt = input_details.init_alt
    init_mach = input_details.init_mach
    mach_final = Ps_init_details.max_mach
    alt_final = Ps_init_details.max_mach_alt
    
    # find energy height of final mach and altitude
    E_w_final = specific_energy(vehicle, alt_final, mach_final)

    # at constant altitude from trajectory starting point, find point of maximum specific excess power; this is the second point of the trajectory
    n1 = 40
    mach_range_temp = np.linspace(init_mach, mach_final, n1)
    Ps_max_init = specific_excess_power_single_alt(vehicle, input_details, init_alt, mach_range_temp)
    i, j = np.where(Ps_max_init.P_s == np.max(Ps_max_init.P_s))
    E_w_init = specific_energy(vehicle, init_alt, init_mach)

    # initialize trajectory data structure. All trajectory points (alt, mach) are stored in this class.
    trajectory = Data()
    trajectory.mach = []
    trajectory.altitude = []
    trajectory.E_w = []

    trajectory.mach.append(init_mach)
    trajectory.mach.append(np.float(mach_range_temp[j]))
    trajectory.altitude.append(init_alt/1000)
    trajectory.altitude.append(init_alt/1000)
    trajectory.E_w.append(np.float(E_w_init))

    # define planet and atmospheric models for gravity and air property calculations
    planet = SUAVE.Attributes.Planets.Earth()
    atmo = SUAVE.Analyses.Atmospheric.US_Standard_1976()

    # increments
    alt_increment = 500 # units [m]
    E_w_increment = 250

    # energy height of second trajectory point (begin climb trajectory from here)
    E_w_current = specific_energy(vehicle, init_alt, trajectory.mach[1])
    trajectory.E_w.append(np.float(E_w_current))
    current_mach = trajectory.mach[-1]
    alt = init_alt

    # initialize toggle. Toggle is used for identifying certain portions of climb trajectory.
    E_w_set = 0

    while E_w_set != 2:

        # if aircraft specific energy is equal to that at final altitude and mach
        if E_w_set == 1:
            E_w_set = 2

        # if specific energy increment goes past that at final altitude and mach, set next energy height to final energy height
        if (E_w_final - E_w_current) < E_w_increment and E_w_set != 2:
            E_w_new = E_w_final
            E_w_set = 1

        # increase energy height by increment
        elif E_w_set == 0:
            E_w_new = E_w_current + E_w_increment
            

        if E_w_set != 2:
            # calculate contour line at new energy height
            n2 = 20
            alt_array = np.linspace((alt - (alt_increment / 2)), (alt + alt_increment + 500), n2)
            mach_array = np.zeros(n2)

            # for given energy height, calculate subsequent speed at a range of altitudes along the contour
            for r, height in enumerate(alt_array):
                gravity_new = planet.compute_gravity(height)
                V = np.sqrt((E_w_new - height) * 2 * gravity_new)
                air_properties_new = atmo.compute_values(height, 0.)
                mach_array[r] = V / air_properties_new.speed_of_sound
            
            # cubic polynomial approximation
            def objective(x, a, b, c, d):
                return a + (b * x) + (c * x**2) + (d * x**3)

            # using cubic polynomial approximation of energy height contour
            popt1, _ = curve_fit(objective, mach_array, alt_array)
            a, b, c, d = popt1
            
            # initialize specific excess power array
            P_s_array = np.zeros(len(mach_array))

            # calculate specific excess power along constant energy height contour
            for o, mach in enumerate(mach_array):
                alt_temp = alt_array[o]
                P_s_array[o] = specific_excess_power_single_point(vehicle, input_details, alt_temp, mach).P_s

            # cubic polynomial approximation of specific excess power curve
            popt2, _ = curve_fit(objective, mach_array, P_s_array)
            e, f, g, h = popt2

            mach_high = mach_array[0]
            mach_low = mach_array[-1]
            mach_delta = 10

            # find point of zero slope along specific excess power curve. Using the bisection method, the mach at the point of zero slope is determined
            while mach_delta > 1e-8:
                P_s_der_high = f + ((2 * g) * mach_high) + ((h * 3) * mach_high**2)

                mach_mid = (mach_high + mach_low) / 2

                P_s_der_mid = f + ((2 * g) * mach_mid) + ((h * 3) * mach_mid**2)

                if P_s_der_mid * P_s_der_high > 0:
                    mach_high = mach_mid
                elif P_s_der_mid * P_s_der_high < 0:
                    mach_low = mach_mid
                
                mach_delta = mach_high - mach_low

            
            current_mach = mach_mid
            
            # calculate altitude along energy height contour for given mach
            alt = objective(current_mach, a, b, c, d)
            E_w_current = specific_energy(vehicle, alt, current_mach)
            
            trajectory.mach.append(current_mach)
            trajectory.altitude.append(alt / 1000)
            trajectory.E_w.append(np.float(E_w_current))
        
        # if final energy height is achieved, follow energy height contour to final mach and altitude
        else:
            n2 = 20
            alt_array = np.linspace(alt, alt_final, n2)
            mach_array = np.zeros(n2)

            for p, height in enumerate(alt_array):
                    gravity_new = planet.compute_gravity(height)
                    V = np.sqrt((E_w_new - height) * 2 * gravity_new)
                    air_properties_new = atmo.compute_values(height, 0.)
                    mach_array[p] = V / air_properties_new.speed_of_sound
                    trajectory.mach.append(mach_array[p])
                    trajectory.altitude.append(height / 1000)
                    trajectory.E_w.append(np.float(E_w_new))
                    E_w_current = E_w_final


    trajectory.mach_final = mach_final
    trajectory.alt_final = alt_final
    
    return trajectory