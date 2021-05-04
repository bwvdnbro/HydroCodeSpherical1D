# set up for the Sedov blastwave test

import get_cmake_command

sedov_options = {
    "rmin_in_au": 0.0,
    "rmax_in_au": 2.64e7,
    "ncell": 1000,
    "gamma": 5.0 / 3.0,
    "maxtime_in_yr": 1.0e5,
    "number_of_snaps": 10,
    "ic": "IC_SEDOV",
    "eos": "EOS_IDEAL",
    "boundaries": "BOUNDARIES_REFLECTIVE",
    "potential": "POTENTIAL_NONE",
    "unit_mass_in_si": 1.99e30,
    "unit_length_in_si": 3.086e19,
    "g_internal": 6.67408e-11,
    "courant_factor": 0.01,
    "riemannsolver_type": "RIEMANNSOLVER_TYPE_EXACT",
}

print((get_cmake_command.get_cmake_command(sedov_options)))
