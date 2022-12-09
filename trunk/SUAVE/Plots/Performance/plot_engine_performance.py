## @ingroup Plots-Performance
# plot_print_engine_performance.py

# Created:  Dec 2022, Sebastian Leblic
# Modified: 
#           

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
import SUAVE
from SUAVE.Core import Units, Data
from SUAVE.Methods.Propulsion import turbofan_sizing, \
    turbojet_small_sizing, \
    turbojet_sizing

import time  # importing library
import datetime  # importing library
import copy


# ----------------------------------------------------------------------
#  Plot max available thrust, 
# ----------------------------------------------------------------------
## @ingroup Plots-Performance
def plot_engine_performance(vehicle, filename='engine_performance.dat', units="SI"):
    """This creates plots showing engine performance through various altitude and speed points.

    Assumptions:
    One network (can be multiple engines). Any engine of the following list can be used:
    - 'Turbofan'
    - 'Turbojet'
    - 'Turbojet_Small'

    Source:
    N/A

    Inputs:
    vehicle.
      tag
      turbofan(), or              Function to compute thrust and fuel burn rate
      tubrojet(), or
      turbojet_small()
      networks.turbofan(or turbojet, turbojet_small).
        design_thrust         [N]
        engine_length         [m]
        thrust.bypass_ratio   [-] 
    filename (optional)       <string> Determines the name of the saved file
    units (optional)          <string> Determines the type of units used in the output, options are imperial and si

    Outputs:
        engine_performance.paramaters.
            max_thrust_available    [N]
            TSFC                    [g/kN.s]
            fuel_flow_rate          [g/s]

    filename                  Saved file with name as above

    Properties Used:
    N/A
    """
    network = vehicle.networks

    for network_type, value in network.items() :

        engine = getattr(vehicle.networks, network_type)
        tag = engine.tag

    imperial = False
    SI = False

    if units.lower() == "imperial":
        imperial = True
    elif units.lower() == "si":
        SI = True
    else:
        print("Incorrect system of units selected - choose 'imperial' or 'SI'")
        return

    alt_min = 0
    alt_max = 14000
    
    mach_min = 0.01
    mach_max = 1.1

    n = 11

    d_isa = 0

    if imperial:
        alt_vec = np.linspace(alt_min, alt_max * 3.28084, n) * Units.ft

    elif SI:
        alt_vec = np.linspace(alt_min, alt_max, n) * Units.m

    mach_vec = np.linspace(mach_min, mach_max, n)
    thrust = np.zeros_like(mach_vec)
    mdot = np.zeros_like(mach_vec)

    # Determining vehicle number of engines
    engine_number = 0.
    
    for network in vehicle.networks:  # may have than one network
        engine_number += network.number_of_engines
    if engine_number == 0:
        raise ValueError("No engine found in the vehicle") 

    # Considering planet and atmosphere of 1st mission segment
    sea_level_gravity = SUAVE.Attributes.Planets.Earth().sea_level_gravity
    atmo = SUAVE.Analyses.Atmospheric.US_Standard_1976()  # mission.segments[0].atmosphere

    atmo_values     = atmo.compute_values(0.,0.)
    
    p0   = atmo_values.pressure
    T0   = atmo_values.temperature
    rho0 = atmo_values.density
    a0   = atmo_values.speed_of_sound
    mu0  = atmo_values.dynamic_viscosity    

    state = Data()
    state.conditions = SUAVE.Analyses.Mission.Segments.Conditions.Aerodynamics()
    state.numerics = SUAVE.Analyses.Mission.Segments.Conditions.Numerics()

    thrust_results = Data()
    thrust_results.thrust = []
    thrust_results.TSFC = []
    thrust_results.fuel_flow_rate = []

    for altitude in alt_vec:
        for idx, mach in enumerate(mach_vec):
            # Computing atmospheric conditions
            atmo_values     = atmo.compute_values(altitude,0)
            
            p   = atmo_values.pressure
            T   = atmo_values.temperature
            rho = atmo_values.density
            a   = atmo_values.speed_of_sound
            mu  = atmo_values.dynamic_viscosity                
            
            T_delta_ISA = T + d_isa
            sigma_disa = (p / p0) / (T_delta_ISA / T0)
            rho_delta_ISA = rho0 * sigma_disa
            a_delta_ISA = atmo.fluid_properties.compute_speed_of_sound(T_delta_ISA)
            speed = mach * a_delta_ISA

            # Getting engine thrust
            state.conditions.freestream.dynamic_pressure = np.array(
                np.atleast_1d(0.5 * rho_delta_ISA * speed * speed))
            state.conditions.freestream.gravity = np.array(np.atleast_1d(sea_level_gravity))
            state.conditions.freestream.velocity = np.array(np.atleast_1d(speed))
            state.conditions.freestream.mach_number = np.array(np.atleast_1d(speed / a_delta_ISA))
            state.conditions.freestream.temperature = np.array(np.atleast_1d(T_delta_ISA))
            state.conditions.freestream.pressure = np.array(np.atleast_1d(p))
            state.conditions.freestream.speed_of_sound = np.array(np.atleast_1d(a_delta_ISA))
            state.conditions.freestream.density = np.array(np.atleast_1d(rho_delta_ISA))
            state.conditions.freestream.dynamic_viscosity = np.array(np.atleast_1d(mu))
            state.conditions.freestream.altitude = np.array(np.atleast_1d(altitude))
            state.conditions.propulsion.throttle = np.array(np.atleast_1d(1.))
            
            if tag =='turbojet_small':
                results = vehicle.networks.turbojet_small.evaluate_thrust(state)  # total thrust
                thrust[idx] = results.thrust_force_vector[0, 0]
                mdot[idx] = results.vehicle_mass_rate[0, 0]
            elif tag =='turbojet':
                results = vehicle.networks.turbojet.evaluate_thrust(state)  # total thrust
                thrust[idx] = results.thrust_force_vector[0, 0]
                mdot[idx] = results.vehicle_mass_rate[0, 0]
            elif tag =='turbofan':
                results = vehicle.networks.turbofan.evaluate_thrust(state)  # total thrust
                thrust[idx] = results.thrust_force_vector[0, 0]
                mdot[idx] = results.vehicle_mass_rate[0, 0]          
            

        if imperial:
            scf = 3600. * mdot / 0.1019715 / thrust
            thrust = np.divide(thrust, engine_number) / Units.lbf
        elif SI:
            scf = mdot / thrust * 1e6
            thrust = np.divide(thrust, engine_number) / Units['N']

        thrust_results.thrust.append(thrust.copy())
        thrust_results.TSFC.append(scf.copy())
        thrust_results.fuel_flow_rate.append(mdot.copy())

    # Plot results
    Mach_range, Alt_range = np.meshgrid(mach_vec, alt_vec)
    
    fig = plt.figure()
    mappable = plt.cm.ScalarMappable(cmap=plt.cm.viridis)
    mappable.set_array(thrust_results.thrust)
    cbar = plt.colorbar(mappable)
    if imperial:
        plt_handle = plt.contourf(Mach_range, Alt_range, thrust_results.thrust, levels=15)
        plt.ylabel('Altitude (ft)')
        cbar.ax.set_ylabel('Thrust (lbf)')   
    
    elif SI:
        plt_handle = plt.contourf(Mach_range, Alt_range/1000, thrust_results.thrust, levels=15)
        plt.ylabel('Altitude (km)')
        cbar.ax.set_ylabel('Thrust (N)')
    plt.xlabel('Mach')
    plt.title('Maximum Thrust Available')

    #plt.savefig('specific_excess_power.pdf')
    #plt.show()

    fig = plt.figure()
    mappable = plt.cm.ScalarMappable(cmap=plt.cm.viridis)
    mappable.set_array(thrust_results.TSFC)
    cbar = plt.colorbar(mappable)
    if imperial:
        plt_handle = plt.contourf(Mach_range, Alt_range, thrust_results.TSFC, levels=15)
        plt.ylabel('Altitude (ft)')
        cbar.ax.set_ylabel('TSFC (lb/lbf.hr)')   
    
    elif SI:
        plt_handle = plt.contourf(Mach_range, Alt_range/1000, thrust_results.TSFC, levels=15)
        plt.ylabel('Altitude (km)')
        cbar.ax.set_ylabel('TSFC (g/kN.s)')

    plt.xlabel('Mach')
    plt.title('Thrust Specific Fuel Consumption (TSFC)')

    #plt.savefig('TSFC.pdf')
    #plt.show()

    fig = plt.figure()
    mappable = plt.cm.ScalarMappable(cmap=plt.cm.viridis)
    mappable.set_array(thrust_results.fuel_flow_rate)
    cbar = plt.colorbar(mappable)
    if imperial:
        plt_handle = plt.contourf(Mach_range, Alt_range, thrust_results.fuel_flow_rate, levels=15)
        plt.ylabel('Altitude (ft)')
        cbar.ax.set_ylabel('Fuel Flow Rate (lb/hr)')   
    
    elif SI:
        plt_handle = plt.contourf(Mach_range, Alt_range/1000, thrust_results.fuel_flow_rate, levels=15)
        plt.ylabel('Altitude (km)')
        cbar.ax.set_ylabel('Fuel Flow Rate (kg/s)')
    
    plt.xlabel('Mach')
    plt.title('Fuel Flow Rate')

    #plt.savefig('fuel_flow_rate.pdf')
    plt.show()

    return


# ----------------------------------------------------------------------
#   Module Test
# ----------------------------------------------------------------------
if __name__ == '__main__':
    print(' Error: No test defined ! ')
