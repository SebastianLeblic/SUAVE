# specific_energy.py
#
# Created: Jan 2023, S. Leblic (AEROCore Labs)
#
# Modified:
#
#
#
#
#

"""
calculates specific energy (or energy height) given a mass, altitude, and velocity

"""

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------

import SUAVE


def specific_energy(vehicle, alt, mach):

    planet = SUAVE.Attributes.Planets.Earth()
    atmo = SUAVE.Analyses.Atmospheric.US_Standard_1976()
    air_properties = atmo.compute_values(alt, 0.)
    gravity = planet.compute_gravity(alt)

    speed = mach * air_properties.speed_of_sound
    mass = vehicle.mass_properties.takeoff
    PE = mass * gravity * alt
    KE = 0.5 * mass * (speed**2)
    E_w = (KE + PE) / (mass * gravity)

    return E_w