# aircraft_energy_min_fuel_climb_template.py
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

Setup file to find energy optimal minimum fuel-to-climb trajectory for given aricraft.

Plots specific excess power map with overlaid trajectory.

"""

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------

import SUAVE
from SUAVE.AEROCore_Additional.energy_climb_trajectory_const_mass import minimum_fuel_to_climb_const_mass
from SUAVE.AEROCore_Additional.fuel_specific_excess_power import fuel_specific_excess_power_range
from SUAVE.AEROCore_Additional.plotting import *
from SUAVE.AEROCore_Additional.aero_forces import construct_aero_map

from SUAVE.Core import Data

import numpy as np
import time

import sys, os
from os import path

# changing directory path to extract vehicle files. Uncomment print basepath statement to see current path and adjust append as necessary
basepath = path.dirname(__file__)
#print("\n basepath:", basepath)
filepath = path.abspath(path.join(basepath, "..", "Vehicles"))
sys.path.append(filepath)
from VEHICLE_NAME_HERE import vehicle_setup


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
    input_details.aero_file = basepath + 'PATH_TO_AERO_COEFFS_FILE_HERE.txt'
    input_details.aero_plotting = False
    input_details.aero_map = construct_aero_map(input_details)
    input_details.d_isa = 0

    # get fuel specific excess power map details
    FPs_details = fuel_specific_excess_power_range(vehicle, input_details, alt_range, mach_range)
    
    # additional initialization for trajectory profiling
    input_details.init_alt = 500 # units [m]
    input_details.init_mach = 0.1

    trajectory_results = minimum_fuel_to_climb_const_mass(vehicle, input_details)
    
    et = time.time()
    elapsed_time = et - st
    print("\n Execution time: ", elapsed_time, 's')
    

    #---------------------------------------------------------------------------------------------
    # PLOTTING
    #---------------------------------------------------------------------------------------------
    FP_s = FPs_details.FP_s
    E_w = FPs_details.E_w
    mach_final = trajectory_results.mach_final
    alt_final = trajectory_results.alt_final
    mass = vehicle.mass_properties.takeoff
    max_alt_plot = max(alt_range) / 1000

    wing_area = vehicle.wings.main_wing.areas.reference
    sealevel_static_thrust = vehicle.networks.turbojet_small.sealevel_static_thrust

    Mach_range, Alt_range = np.meshgrid(mach_range, alt_range)
    plt.rcParams["figure.figsize"] = (19, 10)
    plt.figure('Minimum Fuel to Climb')
    plt_handle = plt.contour(
        Mach_range, Alt_range/1000, E_w, levels=5, colors='k')

    fmt_custom = {}
    for l in plt_handle.levels:
        fmt_custom[l] = "   E/W = " + str(round(l)) + " m  "

    plt.clabel(plt_handle, inline=True, fmt=fmt_custom, fontsize=7)

    Mach_range, Alt_range = np.meshgrid(mach_range, alt_range)
    mappable = plt.cm.ScalarMappable(cmap=plt.cm.viridis)
    mappable.set_array(FP_s)
    plt_handle = plt.contourf(
        Mach_range, Alt_range/1000, FP_s, levels=15)
    cbar = plt.colorbar(mappable)
    
    plt.xlabel('Mach')
    plt.ylabel('Altitude [km]')
    t = 'Vehicle: ' + str(vehicle.tag) + '\nScale: ' + str(round(vehicle.scale, 2)) + '\nWing_area: ' + str(round(wing_area, 3)) + 'm$^{2}$' + '\nAircraft Mass: ' + str(round(mass, 2)) + ' kg' + '\nSea level static thrust: ' + str(round(sealevel_static_thrust, 2)) + ' N'
    plt.text(0.025, 0.8 * max_alt_plot, t, bbox=dict(facecolor='white', alpha=1.))
    cbar.ax.set_ylabel('Fuel specific excess power [m$^{2}$/s$^{2}$]')
    plt.title('Minimum Fuel-to-Climb Trajectory')
    plt.plot(trajectory_results.mach, trajectory_results.altitude, 'w--')
    plt.plot(mach_final, alt_final / 1000, 'ro')

    # plot point of maximum speed
    font = {'color':  'white',
        'weight': 'normal',
        'size': 10,
        }
    t2 = 'Max Mach: ' + str(round(mach_final,3)) + '\nAltitude: ' + str(round(alt_final)) + ' m'
    plt.text(mach_final + 0.01, alt_final / 1000, t2, font)
    plt.ylim([0, max_alt_plot])
    plt.xlim([0, max(mach_range)])
    plt.savefig('Minimum fuel to climb trajectory.pdf')
    plt.show()



    return

def full_setup():

    vehicle = vehicle_setup()

    return vehicle

if __name__ == '__main__': 
    main()



