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

# Configurable parameters

if(NOT DEFINED ic)
  message(STATUS "Using default ic")
  set(ic "IC_FILE")
endif(NOT DEFINED ic)
message(STATUS "IC: ${ic}")

if(NOT DEFINED gamma)
  message(STATUS "Using default gamma")
  set(gamma 1.001)
endif(NOT DEFINED gamma)
message(STATUS "Gamma: ${gamma}")

if(NOT DEFINED dt)
  message(STATUS "Using default dt")
  set(dt 0.00001)
endif(NOT DEFINED dt)
message(STATUS "Dt: ${dt}")

if(NOT DEFINED nstep)
  message(STATUS "Using default nstep")
  set(nstep 200000)
endif(NOT DEFINED nstep)
message(STATUS "Nstep: ${nstep}")

if(NOT DEFINED snapstep)
  message(STATUS "Using default snapstep")
  set(snapstep 200)
endif(NOT DEFINED snapstep)
message(STATUS "Snapstep: ${snapstep}")

if(NOT DEFINED nstep_relax)
  message(STATUS "Using default nstep_relax")
  set(nstep_relax 0)
endif(NOT DEFINED nstep_relax)
message(STATUS "Nstep_relax: ${nstep_relax}")

if(NOT DEFINED eos)
  message(STATUS "Using default eos")
  set(eos "EOS_BONDI")
endif(NOT DEFINED eos)
message(STATUS "EOS: ${eos}")

if(NOT DEFINED bondi_density_exponent)
  message(STATUS "Using default bondi_density_exponent")
  set(bondi_density_exponent 1.2)
endif(NOT DEFINED bondi_density_exponent)
message(STATUS "Bondi_density_exponent: ${bondi_density_exponent}")

if(NOT DEFINED bondi_velocity_exponent)
  message(STATUS "Using default bondi_velocity_exponent")
  set(bondi_velocity_exponent 0.8)
endif(NOT DEFINED bondi_velocity_exponent)
message(STATUS "Bondi_velocity_exponent: ${bondi_velocity_exponent}")

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
