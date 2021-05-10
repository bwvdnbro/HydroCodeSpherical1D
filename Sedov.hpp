/*******************************************************************************
 * This file is part of HydroCodeSpherical1D
 * Copyright (C) 2021 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file Sedov.hpp
 *
 * @brief Set up initial conditions for the Sedov blastwave.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef SEDOV_HPP
#define SEDOV_HPP

#define SEDOV_NCELL_INJECTION 10

/**
 * @brief Initialize the given cells.
 *
 * @param cells Cells to initialize.
 * @param ncell Number of cells.
 */
#if DIMENSIONALITY == DIMENSIONALITY_1D
#define initialize(cells, ncell)                                               \
  _Pragma("omp parallel for") for (unsigned int i = 1; i < ncell + 1; ++i) {   \
    if (i < SEDOV_NCELL_INJECTION + 1) {                                       \
      cells[i]._P = (GAMMA - 1.) * (1. / SEDOV_NCELL_INJECTION) / cells[i]._V; \
    } else {                                                                   \
      cells[i]._P = 1.e-6;                                                     \
    }                                                                          \
    cells[i]._rho = 1.;                                                        \
    cells[i]._u = 0.;                                                          \
    cells[i]._a = 0.;                                                          \
    cells[i]._nfac = 0.;                                                       \
  }
#else
#define initialize(cells, ncell)                                               \
  _Pragma("omp parallel for") for (unsigned int i = 1; i < ncell + 1; ++i) {   \
    if (i < SEDOV_NCELL_INJECTION + 1) {                                       \
      const double Ru = cells[i]._uplim;                                       \
      const double Rl = cells[i]._lowlim;                                      \
      const double Vs = 4. * M_PI / 3. * (Ru * Ru * Ru - Rl * Rl * Rl);        \
      cells[i]._P = (GAMMA - 1.) *                                             \
                    (5.e43 / SEDOV_NCELL_INJECTION / UNIT_ENERGY_IN_SI) / Vs;  \
    } else {                                                                   \
      cells[i]._P = 1.e4 * BOLTZMANN_K_IN_SI / HYDROGEN_MASS_IN_SI /           \
                    UNIT_VELOCITY_IN_SI / UNIT_VELOCITY_IN_SI;                 \
    }                                                                          \
    cells[i]._rho = 4.182e-22 / UNIT_DENSITY_IN_SI;                            \
    cells[i]._u = 0.;                                                          \
    cells[i]._a = 0.;                                                          \
    cells[i]._nfac = 0.;                                                       \
  }
#endif

#endif // SEDOV_HPP
