## @ingroup Analyses-Aerodynamics
# SU2_Euler_Super.py
#
# Created:  Dec 2022, Sebastian Leblic
# Modified: 
#           
#           
#           

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import SUAVE
from SUAVE.Core import Data
from .Markup import Markup
from SUAVE.Analyses import Process
import numpy as np

# The aero methods
from SUAVE.Methods.Aerodynamics.Common import Fidelity_Zero as Common
from SUAVE.Analyses.Aerodynamics.SU2_inviscid_Super_mod import SU2_inviscid_Super_mod

# ----------------------------------------------------------------------
#  Analysis
# ----------------------------------------------------------------------
## @ingroup Analyses-Aerodynamics
class CFD_results(Markup):
    """This uses SU2 to compute lift.

    Assumptions:
    CFD results file exists with. No further aerodynamic analysis is needed.

    Source:
    None
    """     
    def __defaults__(self):
        """This sets the default values and methods for the analysis.

        Assumptions:
        None

        Source:
        N/A

        Inputs:
        None

        Outputs:
        None

        Properties Used:
        N/A
        """              
        self.tag    = 'CFD_results_markup'       
    
        # Correction factors
        settings = self.settings
        settings.trim_drag_correction_factor        = 1.02
        settings.wing_parasite_drag_form_factor     = 1.1
        settings.fuselage_parasite_drag_form_factor = 2.3
        settings.oswald_efficiency_factor           = None
        settings.span_efficiency                    = None
        settings.viscous_lift_dependent_drag_factor = 0.38
        settings.drag_coefficient_increment         = 0.0000
        settings.spoiler_drag_increment             = 0.00 
        settings.maximum_lift_coefficient           = np.inf 
        settings.half_mesh_flag                     = True
        settings.parallel                           = False
        settings.processors                         = 1
        settings.vsp_mesh_growth_ratio              = 1.3
        settings.vsp_mesh_growth_limiting_flag      = False
        settings.recalculate_total_wetted_area      = False
        
        # Build the evaluation process
        compute = self.process.compute
        compute.lift = Process()
        compute.drag = Process()

        # Run SU2
        compute.lift.inviscid                      = SU2_inviscid_Super_mod()
        compute.lift.total                         = Common.Lift.aircraft_total
        
        compute.drag.inviscid                      = SU2_inviscid_Super_mod()
        compute.drag.total                         = Common.Drag.aircraft_total        
        

    def initialize(self):
        """Initializes the surrogate needed for SU2, including building the surface and volume meshes.

        Assumptions:
        Vehicle is available in OpenVSP

        Source:
        N/A

        Inputs:
        None

        Outputs:
        None

        Properties Used:
        None

        """              
        super(CFD_results, self).initialize()
        self.process.compute.lift.inviscid.geometry = self.geometry
        self.process.compute.drag.inviscid.geometry = self.geometry
        
        # Generate the surrogate
        self.process.compute.lift.inviscid.initialize()
        self.process.compute.drag.inviscid.initialize()
        
    finalize = initialize