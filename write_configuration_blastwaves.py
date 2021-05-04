# set up for the Sod shock test

import get_cmake_command

sod_options = {
    "rmin_in_au": 0.0,
    "rmax_in_au": 6.68459e-12,
    "ncell": 1000,
    "gamma": 1.4,
    "maxtime_in_yr": 1.2049721e-9,
    "number_of_snaps": 10,
    "ic": "IC_BLASTWAVES",
    "eos": "EOS_IDEAL",
    "boundaries": "BOUNDARIES_REFLECTIVE",
    "potential": "POTENTIAL_NONE",
    "unit_mass_in_si": 1.0,
    "unit_length_in_si": 1.0,
    "g_internal": 6.67408e-11,
    "courant_factor": 0.01,
    "riemannsolver_type": "RIEMANNSOLVER_TYPE_EXACT",
    "dimensionality": "DIMENSIONALITY_1D",
    "hydro_order": 1,
}

print(get_cmake_command.get_cmake_command(sod_options))
