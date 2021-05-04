# write the configuration command for a starbench run

import get_cmake_command

starbench_options = {
    "rmin_in_au": 0.0,
    "rmax_in_au": 3.0e5,
    "ncell": 1000,
    "gamma": 1.001,
    "maxtime_in_yr": 1.41e5,
    "number_of_snaps": 100,
    "ic": "IC_STARBENCH",
    "eos": "EOS_BONDI",
    "boundaries": "BOUNDARIES_REFLECTIVE",
    "isothermal_temperature_in_k": 100.0,
    "potential": "POTENTIAL_NONE",
    "g_internal": 1.0,
    "mass_point_mass_in_msol": 18.0,
    "bondi_density_in_si": 1.0e-16,
    "bondi_pressure_contrast": 200.0,
    "initial_ionisation_radius_in_au": 6.5e4,
    "unit_mass_in_si": 2.479e31,
    "unit_length_in_si": 1.2e13,
    "ionisation_mode": "IONISATION_MODE_SELF_CONSISTENT",
    "ionisation_transition": "IONISATION_TRANSITION_SMOOTH",
    "ionisation_transition_width_in_au": 3.0e3,
    "courant_factor": 0.05,
    "riemannsolver_type": "RIEMANNSOLVER_TYPE_HLLC",
}

print(get_cmake_command.get_cmake_command(starbench_options))
