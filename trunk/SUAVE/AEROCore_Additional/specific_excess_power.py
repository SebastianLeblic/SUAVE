# specific_excess_power.py
# 
# Created:  Jan 2023, Sebastian Leblic
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

from SUAVE.AEROCore_Additional.aero_forces import aero_forces
from SUAVE.AEROCore_Additional.aero_forces import construct_aero_map
from SUAVE.AEROCore_Additional.calc_propulsion_required import calc_propulsion_required


def specific_excess_power(vehicle, input_details, alt_range, mach_range):

    # unpack
    plotting = input_details.Ps_plotting
    alt_vec = alt_range
    mach_vec = mach_range

    # initializing
    Ps_details = Data()
    Ps_details.P_s = []
    Ps_details.aoa = []
    Ps_details.thrust = []
    Ps_details.TSFC = []
    Ps_details.fuel_flow_rate = []
    Ps_details.E_w = []

    # constructing aero map
    input_details.aero_map = construct_aero_map(input_details)

    # P_s is specific power
    P_s = np.zeros(len(mach_vec))
    aoa = np.zeros(len(mach_vec))
    thrust = np.zeros(len(mach_vec))
    TSFC = np.zeros(len(mach_vec))
    fuel_flow_rate = np.zeros(len(mach_vec))
    # E_w is specific energy
    E_w = np.zeros(len(mach_vec))
    mass = 0
    speed = 0
    gravity = 0
    max_mach = 0

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
            aero_force_results = aero_forces(input_details, vehicle)

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

            
                P_s[j] = ((max_thrust_available - aero_force_results.drag) * speed) / (mass * gravity)
                aoa[j] = aero_force_results.aoa
                thrust[j] = propulsion_results.thrust
                TSFC[j] = propulsion_results.TSFC
                fuel_flow_rate[j] = propulsion_results.fuel_flow_rate

                if P_s[j] < 0:
                    P_s[j] = 0
                    aoa[j] = 0
                    thrust[j] = 0
                    TSFC[j] = 0
                    fuel_flow_rate[j] = 0

                # finding point of maximimum mach. This should be where specific power = 0 at a given altitude.
                if P_s[j] < 50 and P_s[j] >= 0 and mach > 0.8:
                    if mach > max_mach:
                        max_mach = mach
                        max_mach_alt = alt

            else:

                P_s[j] = 0
                aoa[j] = 0
                thrust[j] = 0
                TSFC[j] = 0
                fuel_flow_rate[j] = 0
                speed = aero_force_results.speed

            # Specific energy is calculated
            PE = mass * gravity * alt
            KE = 0.5 * mass * (speed**2)
            E_w[j] = (PE + KE) / (mass * gravity)


        Ps_details.P_s.append(P_s.copy())
        Ps_details.aoa.append(aoa.copy())
        Ps_details.thrust.append(thrust.copy())
        Ps_details.TSFC.append(TSFC.copy())
        Ps_details.fuel_flow_rate.append(fuel_flow_rate.copy())
        Ps_details.E_w.append(E_w.copy())

    Ps_details.max_mach = max_mach
    Ps_details.max_mach_alt = max_mach_alt

    if plotting:
        smoothing = input_details.Ps_plot_smoothing
        save_ps_plots = input_details.save_ps_plots
        plotting_Ps(vehicle, alt_vec, mach_vec, Ps_details, save_ps_plots, smoothing)
    
    return Ps_details
        

def plotting_Ps(vehicle, alt_vec, mach_vec, Ps_details, save_plots, smoothing):

    # unpack
    if smoothing:
        sigma = 0.7
        P_s = gaussian_filter(Ps_details.P_s, sigma)
        aoa = gaussian_filter(Ps_details.aoa, sigma)
        thrust = gaussian_filter(Ps_details.thrust, sigma)
        TSFC = gaussian_filter(Ps_details.TSFC, sigma)
        fuel_flow_rate = gaussian_filter(Ps_details.fuel_flow_rate, sigma)
    
    else:
        P_s = Ps_details.P_s
        aoa = Ps_details.aoa
        thrust = Ps_details.thrust
        TSFC = Ps_details.TSFC
        fuel_flow_rate = Ps_details.fuel_flow_rate
    
    E_w = Ps_details.E_w
    max_mach = Ps_details.max_mach
    max_mach_alt = Ps_details.max_mach_alt
    mass = vehicle.mass_properties.takeoff
    max_alt_plot = max(alt_vec) / 1000

    wing_area = vehicle.wings.main_wing.areas.reference
    sealevel_static_thrust = vehicle.networks.turbojet_small.sealevel_static_thrust

    Mach_range, Alt_range = np.meshgrid(mach_vec, alt_vec)
    plt.rcParams["figure.figsize"] = (19, 10)
    plt.figure('Specific Excess Power')
    plt_handle = plt.contour(
        Mach_range, Alt_range/1000, E_w, levels=5, colors='k')

    fmt_custom = {}
    for l in plt_handle.levels:
        fmt_custom[l] = "   E/W = " + str(round(l)) + "   "

    #plt.ylabel('Altitude (km)')
    #plt.xlabel('Mach')
    plt.clabel(plt_handle, inline=True, fmt=fmt_custom, fontsize=7)

    Mach_range, Alt_range = np.meshgrid(mach_vec, alt_vec)
    mappable = plt.cm.ScalarMappable(cmap=plt.cm.viridis)
    mappable.set_array(P_s)
    plt_handle = plt.contourf(
        Mach_range, Alt_range/1000, P_s, levels=15)
    cbar = plt.colorbar(mappable)
    
    plt.xlabel('Mach')
    plt.ylabel('Altitude [km]')
    t = 'Vehicle: MUFASA_RevA' + '\nScale: 1.0' + '\nWing_area: ' + str(round(wing_area, 3)) + 'm$^{2}$' + '\nAircraft Mass: ' + str(round(mass, 2)) + ' kg' + '\nSea level static thrust: ' + str(round(sealevel_static_thrust, 2)) + ' N'
    plt.text(0.025, 0.8 * max_alt_plot, t, bbox=dict(facecolor='white', alpha=1.))
    cbar.ax.set_ylabel('Specific excess power [m/s]')
    plt.title('Specific Excess Power with Specific Energy Lines')
    plt.plot(max_mach, max_mach_alt / 1000, 'ro')

    # plot point of maximum speed
    font = {'color':  'white',
        'weight': 'normal',
        'size': 10,
        }
    t2 = 'Max Mach: ' + str(round(max_mach,3)) + '\nAltitude: ' + str(round(max_mach_alt)) + ' m'
    plt.text(max_mach + 0.01, max_mach_alt / 1000, t2, font)
    plt.ylim([0, max_alt_plot])
    plt.xlim([0, max(mach_vec)])
    if save_plots:
        plt.savefig('Specific Excess Power.pdf')
    plt.show()

    plt.figure('AoA')
    mappable = plt.cm.ScalarMappable(cmap=plt.cm.viridis)
    mappable.set_array(aoa)
    plt_handle = plt.contourf(
        Mach_range, Alt_range/1000, aoa, levels=15)
    cbar = plt.colorbar(mappable)
    plt.xlabel('Mach')
    plt.ylabel('Altitude (km)')
    cbar.ax.set_ylabel('AoA (deg)')
    plt.title('AoA')
    plt.ylim([0, max_alt_plot])
    plt.xlim([0, max(mach_vec)])
    if save_plots:
        plt.savefig('AoA.pdf')
    plt.show()

    plt.figure('Thrust')
    mappable = plt.cm.ScalarMappable(cmap=plt.cm.viridis)
    mappable.set_array(thrust)
    plt_handle = plt.contourf(
        Mach_range, Alt_range/1000, thrust, levels=15)
    cbar = plt.colorbar(mappable)
    plt.xlabel('Mach')
    plt.ylabel('Altitude (km)')
    cbar.ax.set_ylabel('Thrust [N]')
    plt.title('Thrust')
    plt.ylim([0, max_alt_plot])
    plt.xlim([0, max(mach_vec)])
    if save_plots:
        plt.savefig('Thrust.pdf')
    plt.show()

    plt.figure('TSFC')
    mappable = plt.cm.ScalarMappable(cmap=plt.cm.viridis)
    mappable.set_array(TSFC)
    plt_handle = plt.contourf(
        Mach_range, Alt_range/1000, TSFC, levels=15)
    cbar = plt.colorbar(mappable)
    plt.xlabel('Mach')
    plt.ylabel('Altitude (km)')
    cbar.ax.set_ylabel('TSFC [g/kN.s]')
    plt.title('TSFC')
    plt.ylim([0, max_alt_plot])
    plt.xlim([0, max(mach_vec)])
    if save_plots:
        plt.savefig('TSFC.pdf')
    plt.show()


    plt.figure('Fuel_flow_rate')
    mappable = plt.cm.ScalarMappable(cmap=plt.cm.viridis)
    mappable.set_array(fuel_flow_rate)
    plt_handle = plt.contourf(
        Mach_range, Alt_range/1000, fuel_flow_rate, levels=15)
    cbar = plt.colorbar(mappable)
    plt.xlabel('Mach')
    plt.ylabel('Altitude (km)')
    cbar.ax.set_ylabel('TSFC [g/kN.s]')
    plt.title('Fuel flow rate')
    plt.ylim([0, max_alt_plot])
    plt.xlim([0, max(mach_vec)])
    if save_plots:
        plt.savefig('Fuel_flow_rate.pdf')
    plt.show()
    return