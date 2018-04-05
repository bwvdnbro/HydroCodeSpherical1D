/*******************************************************************************
 * This file is part of HydroCodeSpherical1D
 * Copyright (C) 2018 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file Starbench.hpp
 *
 * @brief Set up initial conditions for the Starbench test.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef STARBENCH_HPP
#define STARBENCH_HPP

/**
 * @brief Initialize the given cells.
 *
 * @param cells Cells to initialize.
 * @param ncell Number of cells.
 */
#define initialize(cells, ncell)                                               \
  _Pragma("omp parallel for") for (unsigned int i = 1; i < ncell + 1; ++i) {   \
    cells[i]._rho = 5.21e-18 / UNIT_DENSITY_IN_SI;                             \
    cells[i]._u = 0.;                                                          \
    cells[i]._P = 1.;                                                          \
    cells[i]._a = 0.;                                                          \
    cells[i]._nfac = 0.;                                                       \
  }

#endif // STARBENCH_HPP
