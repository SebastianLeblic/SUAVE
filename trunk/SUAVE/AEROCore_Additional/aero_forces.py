# @ingroup SUAVE-AEROCore_additional
# aero_forces.py
#
# Created: June 20, 2022 - S. Leblic
# Last modified: 
#           - Dec, 2022 - S. Leblic
#
#

"""
Using aerodynamics coefficients map, interpolates data to find lift and drag
forces.


To-do items:
    1 - Date: June 29, 2022 - Status: not fixed - will need to be fixed when
        plotting mission from takeoff - speeds above approx. 90 m/s are
        feasible at 10m alt.
        Problem: if speed is too low, the cl becomes higher than the aircraft 
        can achieve.
        Solution: allow for contirbution of increased thrust to the lift force
        so as to reduce the required cl, in other words, reducing the required
        lift force produced aerodynamically.
        
"""
# imports
import SUAVE
from SUAVE.Core import Units, Data
import numpy as np
import pylab as plt
from sklearn import gaussian_process
from scipy.optimize import curve_fit
import copy as copy


def aero_forces(input_details, vehicle):
    #print("\n vehicle:", vehicle.wings.main_wing.areas.reference)
    #unpack
    wing_area = vehicle.wings.main_wing.areas.reference
    d_isa = input_details.d_isa
    aero_map = construct_aero_map(input_details)
    altitude = input_details.altitude
    mass = vehicle.mass_properties.takeoff

    # defaults
    thrust_temp = 0
    delta_thrust = 500
    aero_force_set = 0

    # initialize
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

    density = air_properties.density
    mach = input_details.mach
    speed = input_details.mach * air_properties.speed_of_sound

    gravity = planet.compute_gravity(altitude)

    
    while delta_thrust > 1e-6:
        
        # determine aircraft angle attack
        aoa_test, cl_test, cd_test = determine_alpha(vehicle, air_properties, aero_map, input_details, gravity)

        if aoa_test == 100:
            aero_force_set = 1
            delta_thrust = 0
        
        else:
            Thrust1 = 0
            Thrust2 = 0
            
            Thrust1 = 0.5 * (cd_test * density * (speed**2) * wing_area) / np.cos(aoa_test)
            Thrust2 = ((mass * gravity) - ((cl_test * density * (speed**2) * wing_area) / 2)) / np.sin(aoa_test)
            
            lift = 0.5 * cl_test * wing_area * density * (speed ** 2)
            drag = 0.5 * cd_test * wing_area * density * (speed ** 2)
            
            p   = air_properties.pressure
            T   = air_properties.temperature
            a   = air_properties.speed_of_sound
            mu  = air_properties.dynamic_viscosity        
            
            T_delta_ISA = T + d_isa
            sigma_disa = (p / p0) / (T_delta_ISA / T0)
            rho_delta_ISA = density * sigma_disa
            a_delta_ISA = atmo.fluid_properties.compute_speed_of_sound(T_delta_ISA)
            speed = mach * a_delta_ISA
            # adjust thrust intake speed to angle of incidence with freestream
            speed_adj = np.cos(aoa_test) * speed

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
            
            # for min thrust calc
            state.conditions.propulsion.throttle = np.array(np.atleast_1d(0.001))        
            min_thrust_results = vehicle.networks.turbojet_small.evaluate_thrust(state)
            min_thrust = min_thrust_results.thrust_force_vector[0][0]

            # for max thrust calc
            state.conditions.propulsion.throttle = np.array(np.atleast_1d(1.))
            max_thrust_results = vehicle.networks.turbojet_small.evaluate_thrust(state)
            thrust_available = max_thrust_results.thrust_force_vector[0][0]
            relative_error =  np.absolute(Thrust1 - Thrust2) / Thrust2

            if relative_error < 5e-2:
                thrust_total = (Thrust1 + Thrust2) / 2

                throttle_guess = thrust_total / thrust_available
                
                delta_thrust = np.absolute(thrust_total - thrust_temp)
                
                thrust_temp = thrust_total
                

                if thrust_total > thrust_available:
                    aero_force_set = 2
                    delta_thrust = 0
                elif thrust_total < min_thrust:
                    aero_force_set = 3
                    delta_thrust = 0
                else:
                    aero_force_set = 0
                    lift = 0.5 * cl_test * wing_area * density * (speed ** 2)
                    drag = 0.5 * cd_test * wing_area * density * (speed ** 2)
                    aoa = aoa_test / Units.deg
            else:
                print("\n Thrust1:", Thrust1)
                print("\n Thrust2:", Thrust2)
                print('Force balance did not converge!')
                aero_force_set = 1
                delta_thrust = 0
        
  

    
    if aero_force_set == 0:
        
        mass_lifted = lift/gravity

        # setting initial guess for throttle to determine actual available thrust
            
        aero_results = Data()
        aero_results.throttle = throttle_guess
        aero_results.thrust_total = thrust_total
        aero_results.drag = drag
        aero_results.lift = lift
        aero_results.cl = cl_test
        aero_results.cd = cd_test
        aero_results.aoa = aoa
        aero_results.mass_lifted = mass_lifted
        aero_results.speed = speed
        aero_results.gravity = gravity
        aero_results.cl_set = 0
        aero_results.aero_force_set = aero_force_set
        
    else:
        aero_results = Data()
        aero_results.aero_force_set = aero_force_set
        
    return aero_results


def determine_alpha(vehicle, air_properties, aero_map, input_details, gravity ):
    
    # unpack
    density = air_properties.density
    speed = input_details.mach * air_properties.speed_of_sound
    wing_area = vehicle.wings.main_wing.areas.reference
    mass = vehicle.mass_properties.takeoff
    mach = input_details.mach
    
    CL_sur = aero_map.CL_sur
    CD_sur = aero_map.CD_sur
    mach_points = aero_map.mach_points
    AoA_points = aero_map.AoA_points
    mach_percent = np.zeros(len(CL_sur))
    
    for i, elem in enumerate(mach_points):

        curr_el = mach_points[i - 1]
        next_el = elem

        if (mach > curr_el) and (mach < next_el):
            mach_percent = (mach - curr_el)/(next_el - curr_el)
            cl_i = CL_sur[i-1]
            cl_i_next = CL_sur[i] # [i+1]
            cd_i = CD_sur[i-1]
            cd_i_next = CD_sur[i] # [i+1]

    cl_mach = np.zeros(len(cl_i))
    cd_mach = np.zeros(len(cd_i))

    for j, elem in enumerate(cl_i):
        cl_mach[j] = (mach_percent * (cl_i_next[j] - cl_i[j])) + cl_i[j]
        cd_mach[j] = (mach_percent * (cd_i_next[j] - cd_i[j])) + cd_i[j]
    
    # curve fit function of 4th order to cl and cd data from aero map, defined
    # at current mach
    def objective(x, a, b, c, d, e):
        return (a * x) + (b * x**2) + (c * x**3) + (d * x**4) + e

    popt1, _ = curve_fit(objective, AoA_points, cl_mach)
    popt2, _ = curve_fit(objective, AoA_points, cd_mach)
    a, b, c, d, e = popt1
    f, g, h, i, j = popt2
    
    ## iterate to come to solution for alpha
    # going through range of alpha test points, plugging into equation derived
    # from force balance as a function of alpha. Using the simple bisection
    # method, alpha is determined.        
    delta_alpha_test = 500
    aoa_test_high = 0.85
    aoa_test_low = -0.5
    
    counter = 0

    while delta_alpha_test > 1e-10:
        
        # initial bracket
        # high guess
        cl_test_high = objective(aoa_test_high, a, b, c, d, e)
        cd_test_high = objective(aoa_test_high, f, g, h, i, j)

        #solve_high = ((2 * mass * gravity) / (cd_test_high * density * (speed**2) * wing_area)) - (cl_test_high / cd_test_high) - np.tan(aoa_test_high)
        solve_high = ((wing_area * density * (speed**2)) / (2 * mass)) * ((cd_test_high * np.tan(aoa_test_high)) + cl_test_high) - gravity

        # low guess
        cl_test_low = objective(aoa_test_low, a, b, c, d, e)
        cd_test_low = objective(aoa_test_low, f, g, h, i, j)

        #solve_low = ((2 * mass * gravity) / (cd_test_low * density * (speed**2) * wing_area)) - (cl_test_low / cd_test_low) - np.tan(aoa_test_low)
        solve_low = ((wing_area * density * (speed**2)) / (2 * mass)) * ((cd_test_low * np.tan(aoa_test_low)) + cl_test_low) - gravity

        # solve in between
        xns1 = (aoa_test_high + aoa_test_low) / 2
        cl_test = objective(xns1, a, b, c, d, e)
        cd_test = objective(xns1, f, g, h, i, j)

        #solve_temp = ((2 * mass * gravity) / (cd_test * density * (speed**2) * wing_area)) - (cl_test / cd_test) - np.tan(xns1)
        solve_temp = ((wing_area * density * (speed**2)) / (2 * mass)) * ((cd_test * np.tan(xns1)) + cl_test) - gravity

        if (solve_temp * solve_low) > 0:
            aoa_test_low = xns1
            
        if (solve_temp * solve_low) < 0:
            aoa_test_high = xns1
        
        delta_alpha_test = np.absolute(solve_temp)
        
        counter = counter + 1

        if counter > 100:
            delta_alpha_test = 500
            aoa_test_high = 1.5
            aoa_test_low = 1.5
        
        if counter > 500:
            xns1 = 100
            cl_test = 0
            cd_test = 0
            delta_alpha_test = 0
        
        
    aoa = xns1
    cl = cl_test
    cd = cd_test
    
    return aoa, cl, cd

def construct_aero_map(input_details):

    aero_file = input_details.aero_file

    plotting = input_details.plotting
    #plotting = True
    
    data_array = np.loadtxt(aero_file)
    xy = data_array[:, 0:2]
    CL_data = data_array[:, 2:3]
    CD_data = data_array[:, 3:4]

    coefficients = np.hstack([CL_data, CD_data])

    # set prediction model
    regr_cl_sup = gaussian_process._gpr.GaussianProcessRegressor()
    regr_cl_sub = gaussian_process._gpr.GaussianProcessRegressor()
    cl_surrogate_sup = regr_cl_sup.fit(
        xy[xy[:, 1] >= 1.], CL_data[xy[:, 1] >= 1.])
    cl_surrogate_sub = regr_cl_sub.fit(
        xy[xy[:, 1] <= 1.], CL_data[xy[:, 1] <= 1.])
    regr_cd_sup = gaussian_process._gpr.GaussianProcessRegressor()
    regr_cd_sub = gaussian_process._gpr.GaussianProcessRegressor()
    cd_surrogate_sup = regr_cd_sup.fit(
        xy[xy[:, 1] >= 1.], CD_data[xy[:, 1] >= 1.])
    cd_surrogate_sub = regr_cd_sub.fit(
        xy[xy[:, 1] <= 1.], CD_data[xy[:, 1] <= 1.])

    AoA_points = np.linspace(-10, 30, 50) * Units.deg
    mach_points = np.linspace(.005, 1.7, 50)

    AoA_mesh, mach_mesh = np.meshgrid(AoA_points, mach_points)

    CL_sur = np.zeros(np.shape(AoA_mesh))
    CD_sur = np.zeros(np.shape(AoA_mesh))

    for jj in range(len(AoA_points)):
        for ii in range(len(mach_points)):
            if mach_mesh[ii, jj] >= 1.:
                CL_sur[ii, jj] = cl_surrogate_sup.predict(
                    [np.array([AoA_mesh[ii, jj], mach_mesh[ii, jj]])])
                CD_sur[ii, jj] = cd_surrogate_sup.predict(
                    [np.array([AoA_mesh[ii, jj], mach_mesh[ii, jj]])])

            else:
                CL_sur[ii, jj] = cl_surrogate_sub.predict(
                    [np.array([AoA_mesh[ii, jj], mach_mesh[ii, jj]])])
                CD_sur[ii, jj] = cd_surrogate_sub.predict(
                    [np.array([AoA_mesh[ii, jj], mach_mesh[ii, jj]])])
    
    
    if plotting == True:
        mappable = plt.cm.ScalarMappable(cmap=plt.cm.viridis)
        mappable.set_array(CL_sur)
        fig = plt.figure('Coefficient of Lift Surrogate Plot')
        plt_handle = plt.contourf(
            AoA_mesh/Units.deg, mach_mesh, CL_sur, levels=25)
        #plt.clabel(plt_handle, inline=1, fontsize=10)
        cbar = plt.colorbar(mappable)
        plt.scatter(xy[:, 0]/Units.deg, xy[:, 1])
        plt.xlabel('Angle of Attack (deg)')
        plt.ylabel('Mach Number')
        cbar.ax.set_ylabel('Coefficient of Lift')

        plt.show()

        # Stub for plotting drag if implemented:
        mappable = plt.cm.ScalarMappable(cmap=plt.cm.viridis)
        mappable.set_array(CD_sur)
        fig = plt.figure('Coefficient of Drag Surrogate Plot')
        plt_handle = plt.contourf(
            AoA_mesh/Units.deg, mach_mesh, CD_sur, levels=25)
        plt.scatter(xy[:, 0]/Units.deg, xy[:, 1])
        cbar = plt.colorbar(mappable)
        plt.xlabel('Angle of Attack (deg)')
        plt.ylabel('Mach Number')
        cbar.ax.set_ylabel('Coefficient of Drag')

        plt.show()

    aero_map = Data()
    aero_map.CL_sur = CL_sur
    aero_map.CD_sur = CD_sur
    aero_map.AoA_points = AoA_points
    aero_map.mach_points = mach_points
        
    return aero_map
