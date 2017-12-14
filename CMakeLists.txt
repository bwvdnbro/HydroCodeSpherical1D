################################################################################
# This file is part of HydroCodeSpherical1D
# Copyright (C) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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

cmake_minimum_required(VERSION 2.8)

project(HydroCodeSpherical1D)

set(SOURCES
    main_spherical.cpp

    RiemannSolver.hpp
)

# macro that checks if a configuration value was set on the command line
# if not, the default value is set and a warning is printed
# we also print the actual value that is used
macro(check_configuration_option option_name option_default_value)
  if(NOT DEFINED ${option_name})
    message(STATUS "Using default ${option_name}")
    set(${option_name} ${option_default_value})
  endif(NOT DEFINED ${option_name})
  message(STATUS "${option_name}: ${${option_name}}")
endmacro(check_configuration_option)

# Configurable parameters

check_configuration_option(rmin_in_au 10.)
check_configuration_option(rmax_in_au 100.)
check_configuration_option(ncell 2700)
check_configuration_option(gamma 1.001)
check_configuration_option(maxtime_in_yr 80.)
check_configuration_option(number_of_snaps 4000)
check_configuration_option(ic "IC_FILE")
check_configuration_option(eos "EOS_BONDI")
check_configuration_option(boundaries "BOUNDARIES_BONDI")
check_configuration_option(isothermal_temperature_in_k 500.)
check_configuration_option(potential "POTENTIAL_POINT_MASS")
check_configuration_option(g_internal 1.)
check_configuration_option(mass_point_mass_in_msol 18.)
check_configuration_option(bondi_density_in_si 1.e-16)
check_configuration_option(bondi_pressure_contrast 32.)
check_configuration_option(ic_file_name "ic.dat")
check_configuration_option(initial_ionisation_radius_in_au 30.)
check_configuration_option(unit_mass_in_si 2.479e31)
check_configuration_option(unit_length_in_si 1.2e13)
check_configuration_option(ionisation_mode "IONISATION_MODE_SELF_CONSISTENT")
check_configuration_option(ionisation_transition "IONISATION_TRANSITION_SMOOTH")
check_configuration_option(ionisation_transition_width_in_au 5.)

configure_file(${PROJECT_SOURCE_DIR}/Parameters.hpp.in
               ${PROJECT_SOURCE_DIR}/Parameters.hpp @only)

find_package(OpenMP)
if(OPENMP_FOUND)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
else(OPENMP_FOUND)
  message(ERROR "OpenMP not found!")
endif(OPENMP_FOUND)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ffast-math -ftree-vectorize")

add_executable(HydroCodeSpherical1D ${SOURCES})
