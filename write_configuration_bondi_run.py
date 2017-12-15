# get the configuration command for a Bondi stability test run

import get_cmake_command

bondi_options = {
"rmin_in_au": 10.,
"rmax_in_au": 100.,
"ncell": 2700,
"gamma": 1.001,
"maxtime_in_yr": 80.,
"number_of_snaps": 4000,
"ic": "IC_FILE",
"eos": "EOS_BONDI",
"boundaries": "BOUNDARIES_BONDI",
"isothermal_temperature_in_k": 500.,
"potential": "POTENTIAL_POINT_MASS",
"g_internal": 1.,
"mass_point_mass_in_msol": 18.,
"bondi_density_in_si": 1.e-16,
"bondi_pressure_contrast": 32.,
"ic_file_name": "bondi_ic.dat",
"initial_ionisation_radius_in_au": 30.,
"unit_mass_in_si": 2.479e31,
"unit_length_in_si": 1.2e13,
"ionisation_mode": "IONISATION_MODE_SELF_CONSISTENT",
"ionisation_transition": "IONISATION_TRANSITION_SMOOTH",
"ionisation_transition_width_in_au": 5.,
"courant_factor": 0.05,
"riemannsolver_type": "RIEMANNSOLVER_TYPE_HLLC"
}

print get_cmake_command.get_cmake_command(bondi_options)
