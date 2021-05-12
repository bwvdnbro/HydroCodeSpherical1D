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

/**
 * @file Boundaries.hpp
 *
 * @brief Boundary conditions.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef BOUNDARIES_HPP
#define BOUNDARIES_HPP

#if BOUNDARIES == BOUNDARIES_BONDI

// Bondi specific boundary conditions are in Bondi.hpp
#include "Bondi.hpp"

#else // BOUNDARIES == BOUNDARIES_BONDI

/**
 * @brief Initialize variables used for the boundary conditions.
 *
 * Open or reflective boundaries don't have associated variables, so this method
 * does nothing.
 */
#define boundary_conditions_initialize()

/**
 * @brief Apply boundary conditions after the primitive variable conversion.
 */
#if BOUNDARIES == BOUNDARIES_OPEN
#define boundary_conditions_primitive_variables()                              \
  /* just mirror the values across the boundary */                             \
  cells[0]._dt = cells[1]._dt;                                                 \
  cells[0]._rho = cells[1]._rho;                                               \
  cells[0]._u = cells[1]._u;                                                   \
  cells[0]._P = cells[1]._P;                                                   \
  cells[ncell + 1]._dt = cells[ncell]._dt;                                     \
  cells[ncell + 1]._rho = cells[ncell]._rho;                                   \
  cells[ncell + 1]._u = cells[ncell]._u;                                       \
  cells[ncell + 1]._P = cells[ncell]._P;
#elif BOUNDARIES == BOUNDARIES_REFLECTIVE
#define boundary_conditions_primitive_variables()                              \
  /* mirror the variables and reverse the sign of the velocity */              \
  cells[0]._dt = cells[1]._dt;                                                 \
  cells[0]._rho = cells[1]._rho;                                               \
  cells[0]._u = -cells[1]._u;                                                  \
  cells[0]._P = cells[1]._P;                                                   \
  cells[ncell + 1]._dt = cells[ncell]._dt;                                     \
  cells[ncell + 1]._rho = cells[ncell]._rho;                                   \
  cells[ncell + 1]._u = -cells[ncell]._u;                                      \
  cells[ncell + 1]._P = cells[ncell]._P;
#endif

/**
 * @brief Apply boundary conditions after the gradient computation.
 */
#if BOUNDARIES == BOUNDARIES_OPEN
#define boundary_conditions_gradients()                                        \
  /* just mirror the gradients across the boundary */                          \
  cells[0]._grad_rho = cells[1]._grad_rho;                                     \
  cells[0]._grad_u = cells[1]._grad_u;                                         \
  cells[0]._grad_P = cells[1]._grad_P;                                         \
  cells[ncell + 1]._grad_rho = cells[ncell]._grad_rho;                         \
  cells[ncell + 1]._grad_u = cells[ncell]._grad_u;                             \
  cells[ncell + 1]._grad_P = cells[ncell]._grad_P;
#elif BOUNDARIES == BOUNDARIES_REFLECTIVE
#define boundary_conditions_gradients()                                        \
  /* reverse the sign of the gradients and do some magic to get an accurate    \
     gradient for the velocities */                                            \
  cells[0]._grad_rho = -cells[1]._grad_rho;                                    \
  cells[0]._grad_u =                                                           \
      -cells[1]._grad_u - cells[1]._grad_u -                                   \
      4. * cells[0]._u / (cells[0]._midpoint - cells[1]._midpoint);            \
  cells[0]._grad_P = -cells[1]._grad_P;                                        \
  cells[ncell + 1]._grad_rho = -cells[ncell]._grad_rho;                        \
  cells[ncell + 1]._grad_u =                                                   \
      -cells[ncell]._grad_u -                                                  \
      4. * cells[ncell]._u /                                                   \
          (cells[ncell]._midpoint - cells[ncell + 1]._midpoint);               \
  cells[ncell + 1]._grad_P = -cells[ncell]._grad_P;

#endif

/**
 * @brief Apply boundary conditions right before the flux calculation, in
 * case of the left boundary.
 */
#if BOUNDARIES == BOUNDARIES_REFLECTIVE
#define boundary_left(rhoL, uL, PL, rhoR, uR, PR)                              \
  rhoL = rhoR;                                                                 \
  uL = -uR;                                                                    \
  PL = PR;
#else
#define boundary_left(rhoL, uL, PL, rhoR, uR, PR)                              \
  rhoL = rhoR;                                                                 \
  uL = uR;                                                                     \
  PL = PR;
#endif

/**
 * @brief Apply boundary conditions right before the flux calculation, in
 * case of the right boundary.
 */
#if BOUNDARIES == BOUNDARIES_REFLECTIVE
#define boundary_right(rhoL, uL, PL, rhoR, uR, PR)                             \
  rhoR = rhoL;                                                                 \
  uR = -uL;                                                                    \
  PR = PL;
#else
#define boundary_right(rhoL, uL, PL, rhoR, uR, PR)                             \
  rhoR = rhoL;                                                                 \
  uR = uL;                                                                     \
  PR = PL;
#endif

#endif // BOUNDARIES == BOUNDARIES_BONDI

#endif // BOUNDARIES_HPP
