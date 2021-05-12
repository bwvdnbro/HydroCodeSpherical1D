#! /usr/bin/python

################################################################################
# This file is part of HydroCodeSpherical1D
# Copyright (C) 2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
#
# HydroCodeSpherical1D is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# HydroCodeSpherical1D is distributed in the hope that it will be useful,
# but WITOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with HydroCodeSpherical1D. If not, see <http://www.gnu.org/licenses/>.
################################################################################

##
# @file get_cmake_command.py
#
# @brief Script that generates the cmake command necessary to configure the code
# with a specific configuration.
#
# @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
##

# default options: do not touch this!

configuration_options = {
    "rmin_in_au": 10.0,
    "rmax_in_au": 100.0,
    "ncell": 2700,
    "gamma": 1.001,
    "maxtime_in_yr": 80.0,
    "number_of_snaps": 4000,
    "ic": "IC_FILE",
    "eos": "EOS_BONDI",
    "boundaries": "BOUNDARIES_BONDI",
    "isothermal_temperature_in_k": 500.0,
    "potential": "POTENTIAL_POINT_MASS",
    "g_internal": 1.0,
    "mass_point_mass_in_msol": 18.0,
    "bondi_density_in_si": 1.0e-16,
    "bondi_pressure_contrast": 32.0,
    "ic_file_name": "ic.dat",
    "initial_ionisation_radius_in_au": 30.0,
    "unit_mass_in_si": 2.479e31,
    "unit_length_in_si": 1.2e13,
    "ionisation_mode": "IONISATION_MODE_SELF_CONSISTENT",
    "ionisation_transition": "IONISATION_TRANSITION_SMOOTH",
    "ionisation_transition_width_in_au": 5.0,
    "courant_factor": 0.05,
    "riemannsolver_type": "RIEMANNSOLVER_TYPE_HLLC",
    "dimensionality": "DIMENSIONALITY_3D",
    "hydro_order": 2,
    "status_update_interval": 10.0,
    "cooling": "COOLING_NONE",
    "minimum_temperature_in_k": 0.0,
}

##
# @brief Generate the cmake command to configure the code with a specific
# configuration.
#
# @param custom_options Configuration options that should replace their
# respective default values.
# @param folder Folder where the CMakeLists.txt file is located.
# @return CMake configuration command.
##
def get_cmake_command(custom_options={}, folder=".."):
    global configuration_options

    command = "cmake -DCMAKE_BUILD_TYPE=Release"
    for option in configuration_options:
        if option in custom_options:
            value = custom_options[option]
            custom_options[option] = "read"
        else:
            value = configuration_options[option]
        command += " -D{0}={1}".format(option, value)
    command += " " + folder

    # check that all custom options were actually used
    for option in custom_options:
        if not custom_options[option] == "read":
            print("Unknown option:", option)
            exit()

    return command


if __name__ == "__main__":
    print(get_cmake_command())
