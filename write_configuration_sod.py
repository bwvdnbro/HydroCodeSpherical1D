import get_cmake_command

sod_options = {
"rmin_in_au": 0.,
"rmax_in_au": 6.68459e-12,
"ncell": 1000,
"gamma": 5. / 3.,
"maxtime_in_yr": 6.342e-9,
"number_of_snaps": 10,
"ic": "IC_SOD",
"eos": "EOS_IDEAL",
"boundaries": "BOUNDARIES_OPEN",
"potential": "POTENTIAL_NONE",
"unit_mass_in_si": 1.,
"unit_length_in_si": 1.,
"g_internal": 6.67408e-11,
"courant_factor": 0.01
}

print get_cmake_command.get_cmake_command(sod_options)
