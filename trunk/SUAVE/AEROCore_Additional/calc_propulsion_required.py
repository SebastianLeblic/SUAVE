# @ingroup SUAVE-AEROCore_additional
# calc_propulsion_required.py
#
# Created: June 20, 2022 - S. Leblic
# Last modified: 
#           - Dec, 2022 - S. Leblic
#
#

"""
Takes thrust required and outputs propulsion parameters to meet thrust needed.

Inputs:
    - thrust_required   [N]
    - throttle_guess    [%] (initial throttle position estimate)

Outputs:
    - throttle          [%]
    - fuel_flow_rate    [kg/s]
    - TSFC              [g/kN.s]

"""

# imports
import SUAVE
from SUAVE.Core import Units, Data
import numpy as np
import pylab as plt
import copy as copy

def calc_propulsion_required(thrust_requirements, vehicle, input_details):
    # unpack
    altitude = input_details.altitude
    d_isa = input_details.d_isa
    aoa = thrust_requirements.aoa
    throttle_guess = thrust_requirements.throttle
    thrust_required = thrust_requirements.thrust

    # initialize
    delta_thrust = 500
    sea_level_gravity = SUAVE.Attributes.Planets.Earth().sea_level_gravity
    planet = SUAVE.Attributes.Planets.Earth()
    atmo = SUAVE.Analyses.Atmospheric.US_Standard_1976()
    air_properties = atmo.compute_values(altitude, 0.)
    atmo_values     = atmo.compute_values(0.,0.)
    
    p0   = atmo_values.pressure
    T0   = atmo_values.temperature
    rho0 = atmo_values.density
    a0   = atmo_values.speed_of_sound
    mu0  = atmo_values.dynamic_viscosity 

    mach = input_details.mach
    speed = input_details.mach * air_properties.speed_of_sound

    gravity = planet.compute_gravity(altitude)

    p   = air_properties.pressure
    T   = air_properties.temperature
    a   = air_properties.speed_of_sound
    mu  = air_properties.dynamic_viscosity
    rho = air_properties.density       
    
    T_delta_ISA = T + d_isa
    sigma_disa = (p / p0) / (T_delta_ISA / T0)
    rho_delta_ISA = rho * sigma_disa
    a_delta_ISA = atmo.fluid_properties.compute_speed_of_sound(T_delta_ISA)
    speed = mach * a_delta_ISA
    ref_length = vehicle.fuselages.fuselage.lengths.total
    reynolds_number = rho * speed * ref_length / mu

    # adjust thrust intake speed to angle of incidence with freestream
    speed_adj = np.cos(aoa) * speed

    state = Data()
    state.conditions = SUAVE.Analyses.Mission.Segments.Conditions.Aerodynamics()
    state.numerics = SUAVE.Analyses.Mission.Segments.Conditions.Numerics()
    state.conditions.freestream.dynamic_pressure = np.array(np.atleast_1d(0.5 * rho_delta_ISA * speed_adj * speed_adj))
    state.conditions.freestream.gravity = np.array(np.atleast_1d(gravity))
    state.conditions.freestream.velocity = np.array(np.atleast_1d(speed_adj))
    state.conditions.freestream.mach_number = np.array(np.atleast_1d(speed_adj / a_delta_ISA))
    state.conditions.freestream.temperature = np.array(np.atleast_1d(T_delta_ISA))
    state.conditions.freestream.pressure = np.array(np.atleast_1d(p))
    state.conditions.freestream.speed_of_sound = np.array(np.atleast_1d(a_delta_ISA))
    state.conditions.freestream.density = np.array(np.atleast_1d(rho_delta_ISA))
    state.conditions.freestream.dynamic_viscosity = np.array(np.atleast_1d(mu))
    state.conditions.freestream.altitude = np.array(np.atleast_1d(altitude))
    state.conditions.freestream.reynolds_number = np.array(np.atleast_1d(reynolds_number))

    other = Data()
    other.isentropic_expansion_factor = 1.4
    other.mach = mach
    
    # calc max thrust properties
    state.conditions.propulsion.throttle = np.array(np.atleast_1d(1.))        
    max_thrust_results = vehicle.networks.turbojet_small.evaluate_thrust(state)
    max_thrust = max_thrust_results.thrust_force_vector[0][0]
    max_fuel_flow = max_thrust_results.vehicle_mass_rate

    # find thrust output parameters
    while delta_thrust > 1e-10:

        state.conditions.propulsion.throttle = np.array(np.atleast_1d(throttle_guess))        
        thrust_results = vehicle.networks.turbojet_small.evaluate_thrust(state)
        thrust = thrust_results.thrust_force_vector[0][0]
        epsilon = thrust / max_thrust
        
        # fuel flow rate adjusted to provide more realistic results. (B. Dalman, eq 4.12)
        fuel_flow_rate = max_fuel_flow * 0.071 * np.exp(2.651 * (np.log(epsilon) - np.log(0.012)) / 4.473)

        delta_thrust = np.absolute(thrust - thrust_required) / thrust_required

        if delta_thrust > 1e-10:

            if thrust > thrust_required:
                throttle_temp = throttle_guess
                throttle_guess = throttle_temp * (1 - (np.absolute(thrust - thrust_required)) / thrust_required)
            
            elif thrust < thrust_required:
                throttle_temp = throttle_guess
                throttle_guess = throttle_temp * (1 + (np.absolute(thrust - thrust_required)) / thrust_required)
        
        else:
            propulsion_results = Data()
            propulsion_results.throttle = throttle_guess
            propulsion_results.TSFC = (fuel_flow_rate / 1000) / (thrust / 1000)
            propulsion_results.thrust = thrust
            propulsion_results.fuel_flow_rate = fuel_flow_rate

    propulsion_results.max_thrust = max_thrust

    return propulsion_results