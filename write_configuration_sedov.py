# set up for the Sedov blastwave test

import get_cmake_command
import sys

dim = 3
if len(sys.argv) > 1:
    dim = int(sys.argv[1])

if dim == 3:
    #    sedov_options = {
    #        "rmin_in_au": 0.0,
    #        "rmax_in_au": 6.68459e-12,
    #        "ncell": 1000,
    #        "gamma": 5.0 / 3.0,
    #        "maxtime_in_yr": 6.342e-9,
    #        "number_of_snaps": 10,
    #        "ic": "IC_SEDOV",
    #        "eos": "EOS_IDEAL",
    #        "boundaries": "BOUNDARIES_REFLECTIVE",
    #        "potential": "POTENTIAL_NONE",
    #        "unit_mass_in_si": 1.0,
    #        "unit_length_in_si": 1.0,
    #        "g_internal": 6.67408e-11,
    #        "courant_factor": 0.1,
    #        "riemannsolver_type": "RIEMANNSOLVER_TYPE_HLLC",
    #        "dimensionality": "DIMENSIONALITY_3D",
    #        "hydro_order": 2,
    #    }
    sedov_options = {
        "rmin_in_au": 0.0,
        "rmax_in_au": 2.64e7,
        "ncell": 1000,
        "gamma": 5.0 / 3.0,
        "maxtime_in_yr": 1.0e6,
        "number_of_snaps": 10,
        "ic": "IC_SEDOV",
        "eos": "EOS_IDEAL",
        "boundaries": "BOUNDARIES_REFLECTIVE",
        "potential": "POTENTIAL_NONE",
        "unit_mass_in_si": 1.99e30,
        "unit_length_in_si": 3.086e19,
        "g_internal": 6.67408e-11,
        "courant_factor": 0.1,
        "riemannsolver_type": "RIEMANNSOLVER_TYPE_EXACT",
        "dimensionality": "DIMENSIONALITY_3D",
    }
else:
    sedov_options = {
        "rmin_in_au": 0.0,
        "rmax_in_au": 6.68459e-12,
        "ncell": 1000,
        "gamma": 5.0 / 3.0,
        "maxtime_in_yr": 6.342e-9,
        "number_of_snaps": 10,
        "ic": "IC_SEDOV",
        "eos": "EOS_IDEAL",
        "boundaries": "BOUNDARIES_REFLECTIVE",
        "potential": "POTENTIAL_NONE",
        "unit_mass_in_si": 1.0,
        "unit_length_in_si": 1.0,
        "g_internal": 6.67408e-11,
        "courant_factor": 0.1,
        "riemannsolver_type": "RIEMANNSOLVER_TYPE_HLLC",
        "dimensionality": "DIMENSIONALITY_1D",
        "hydro_order": 2,
    }


print((get_cmake_command.get_cmake_command(sedov_options)))
