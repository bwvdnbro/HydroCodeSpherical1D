/*******************************************************************************
 * This file is part of HydroCodeSpherical1D
 * Copyright (C) 2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
 *
 * HydroCodeSpherical1D is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * HydroCodeSpherical1D is distributed in the hope that it will be useful,
 * but WITOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with HydroCodeSpherical1D. If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/

/*! @brief Minimum radius (in internal units of L). */
#define RMIN 0.1
/*! @brief Maximum radius (in internal units of L). */
#define RMAX 0.475

/*! @brief Number of spherical shells (1D "cells") in between RMIN and RMAX. */
#define NCELL 2000
/*! @brief Adiabatic index. */
#define GAMMA 1.001
/*! @brief Fixed time step (in internal units of T). */
#define DT 0.00001
/*! @brief Number of time steps. */
#define NSTEP 200000
/*! @brief Number of time steps between subsequent snapshot dumps. */
#define SNAPSTEP 200

/*! @brief Number of time steps before switching on ionization. */
#define NSTEP_RELAX 0

/*! @brief Width of the ionization transition region. */
#define IONIZATION_TRANSITION_WIDTH 0.05

/*! @brief Choice of boundary conditions. */
#define BOUNDARIES BOUNDARIES_BONDI

/*! @brief Choice of equation of state. */
#define EOS EOS_BONDI

/*! @brief Isothermal internal energy (if EOS_ISOTHERMAL is selected, in
 *  internal units of L^2 T^-2). */
#define ISOTHERMAL_U 1.145070026

/*! @brief Constant temperature in the neutral region (in K). */
#define ISOTHERMAL_TEMPERATURE 100.

/*! @brief Choice of external potential. */
#define POTENTIAL POTENTIAL_POINT_MASS

/*! @brief Value for Newton's gravity constant (in internal units of L^3 M^-1
 *  T^-2). */
#define G 1.

/*! @brief Mass of the point mass (if POTENTIAL_POINT_MASS is selected, in
 *  internal units of M). */
#define MASS_POINT_MASS 1.5

/*! @brief Switch off gradients (reduce to first order hydro scheme). */
//#define NO_GRADIENTS

/*! @brief Choice of initial conditions. */
#define IC IC_FILE

/*! @brief Approximate initial ionization radius (in internal units of L). */
#define INITIAL_IONIZATION_RADIUS 0.3

/*! @brief Mass unit in kg. */
#define UNIT_MASS_IN_SI 2.479e31

/*! @brief Length unit in m. */
#define UNIT_LENGTH_IN_SI 1.2e13
