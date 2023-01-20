# plotting.py
#
# Created: Jan 2023, S. Leblic (AEROCore Labs)
#
# Modified: 
#
#
#
#

"""
Various plotting functions can called here.

"""
#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage.filters import gaussian_filter

def plotting_Ps(vehicle, alt_vec, mach_vec, Ps_details, save_plots, smoothing):
    
    """
    plots specific excess power envelope, aoa, thrust, TSFC, fuel_flow_Rate
    """
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
    t = 'Vehicle: ' + str(vehicle.tag) + '\nScale: ' + str(round(vehicle.scale, 2)) + '\nWing_area: ' + str(round(wing_area, 3)) + 'm$^{2}$' + '\nAircraft Mass: ' + str(round(mass, 2)) + ' kg' + '\nSea level static thrust: ' + str(round(sealevel_static_thrust, 2)) + ' N'
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
    plt.plot(max_mach, max_mach_alt / 1000, 'ro')
    plt.ylim([0, max_alt_plot])
    plt.xlim([0, max(mach_vec)])
    if save_plots:
        plt.savefig('Fuel_flow_rate.pdf')
    plt.show()

    return

def plotting_FPs(vehicle, alt_vec, mach_vec, FPs_details, save_plots, smoothing):

    # unpack
    if smoothing:
        sigma = 0.7
        FP_s = gaussian_filter(FPs_details.FP_s, sigma)
        aoa = gaussian_filter(FPs_details.aoa, sigma)
        thrust = gaussian_filter(FPs_details.thrust, sigma)
        TSFC = gaussian_filter(FPs_details.TSFC, sigma)
        fuel_flow_rate = gaussian_filter(FPs_details.fuel_flow_rate, sigma)
    
    else:
        FP_s = FPs_details.FP_s
        aoa = FPs_details.aoa
        thrust = FPs_details.thrust
        TSFC = FPs_details.TSFC
        fuel_flow_rate = FPs_details.fuel_flow_rate
    
    E_w = FPs_details.E_w
    mass = vehicle.mass_properties.takeoff
    max_alt_plot = max(alt_vec) / 1000

    wing_area = vehicle.wings.main_wing.areas.reference
    sealevel_static_thrust = vehicle.networks.turbojet_small.sealevel_static_thrust

    Mach_range, Alt_range = np.meshgrid(mach_vec, alt_vec)
    plt.rcParams["figure.figsize"] = (19, 10)
    plt.figure('Fuel Specific Excess Power')
    plt_handle = plt.contour(
        Mach_range, Alt_range/1000, E_w, levels=5, colors='k')

    fmt_custom = {}
    for l in plt_handle.levels:
        fmt_custom[l] = "   E/W = " + str(round(l)) + "   "

    plt.clabel(plt_handle, inline=True, fmt=fmt_custom, fontsize=7)

    Mach_range, Alt_range = np.meshgrid(mach_vec, alt_vec)
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
    plt.title('Fuel Specific Excess Power with Specific Energy Lines')

    plt.ylim([0, max_alt_plot])
    plt.xlim([0, max(mach_vec)])
    if save_plots:
        plt.savefig('Fuel Specific Excess Power.pdf')
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


def plotting_FFPs(vehicle, alt_vec, mach_vec, FFPs_details, save_plots, smoothing):

    # unpack
    if smoothing:
        sigma = 0.7
        FFP_s = gaussian_filter(FFPs_details.FFP_s, sigma)
        aoa = gaussian_filter(FFPs_details.aoa, sigma)
        thrust = gaussian_filter(FFPs_details.thrust, sigma)
        TSFC = gaussian_filter(FFPs_details.TSFC, sigma)
        fuel_flow_rate = gaussian_filter(FFPs_details.fuel_flow_rate, sigma)
    
    else:
        FFP_s = FFPs_details.FFP_s
        aoa = FFPs_details.aoa
        thrust = FFPs_details.thrust
        TSFC = FFPs_details.TSFC
        fuel_flow_rate = FFPs_details.fuel_flow_rate
    
    E_w = FFPs_details.E_w
    mass = vehicle.mass_properties.takeoff
    max_alt_plot = max(alt_vec) / 1000

    wing_area = vehicle.wings.main_wing.areas.reference
    sealevel_static_thrust = vehicle.networks.turbojet_small.sealevel_static_thrust

    Mach_range, Alt_range = np.meshgrid(mach_vec, alt_vec)
    plt.rcParams["figure.figsize"] = (19, 10)
    plt.figure('Fuel Flow Specific Excess Power')
    plt_handle = plt.contour(
        Mach_range, Alt_range/1000, E_w, levels=5, colors='k')

    fmt_custom = {}
    for l in plt_handle.levels:
        fmt_custom[l] = "   E/W = " + str(round(l)) + "   "

    plt.clabel(plt_handle, inline=True, fmt=fmt_custom, fontsize=7)

    Mach_range, Alt_range = np.meshgrid(mach_vec, alt_vec)
    mappable = plt.cm.ScalarMappable(cmap=plt.cm.viridis)
    mappable.set_array(FFP_s)
    plt_handle = plt.contourf(
        Mach_range, Alt_range/1000, FFP_s, levels=15)
    cbar = plt.colorbar(mappable)
    
    plt.xlabel('Mach')
    plt.ylabel('Altitude [km]')
    t = 'Vehicle: ' + str(vehicle.tag) + '\nScale: ' + str(round(vehicle.scale, 2)) + '\nWing_area: ' + str(round(wing_area, 3)) + 'm$^{2}$' + '\nAircraft Mass: ' + str(round(mass, 2)) + ' kg' + '\nSea level static thrust: ' + str(round(sealevel_static_thrust, 2)) + ' N'
    plt.text(0.025, 0.8 * max_alt_plot, t, bbox=dict(facecolor='white', alpha=1.))
    cbar.ax.set_ylabel('Fuel flow specific excess power [m/kg]')
    plt.title('Fuel Flow Specific Excess Power with Specific Energy Lines')

    plt.ylim([0, max_alt_plot])
    plt.xlim([0, max(mach_vec)])
    if save_plots:
        plt.savefig('Fuel Flow Specific Excess Power.pdf')
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