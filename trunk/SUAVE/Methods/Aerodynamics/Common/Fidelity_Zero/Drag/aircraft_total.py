## @ingroup Methods-Aerodynamics-Common-Fidelity_Zero-Drag
# aircraft_total.py
# 
# Created:  Jun 2022, S. Leblic,
# Modified: 
#      
#

# ----------------------------------------------------------------------
#  Aircraft Total
# ----------------------------------------------------------------------

def aircraft_total(state,settings,geometry):
    """Returns total aircraft drag and stores values

    Assumptions:
    None

    Source:
    None

    Inputs:
    state.conditions.aerodynamics.drag_coefficient    [Unitless]

    Outputs:
    aircraft_lift_total (drag coefficient)            [Unitless]

    Properties Used:
    N/A
    """      
    # print(state.conditions.aerodynamics.inviscid_drag_coefficient)
    aircraft_drag_total = state.conditions.aerodynamics.inviscid_drag_coefficient


    return aircraft_drag_total