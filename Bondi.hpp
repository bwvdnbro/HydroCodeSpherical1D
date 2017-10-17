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
 * @file Bondi.hpp
 *
 * @brief Bondi solution expressions.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef BONDI_HPP
#define BONDI_HPP

#include "Cell.hpp"
#include "SafeParameters.hpp"

#include <cmath>

/**
 * @brief Get the value of the Bondi density at the given inverse radius.
 *
 * @param rinv Inverse radius (in internal units of L^-1).
 * @return Value of the density (in internal units of M L^-3).
 */
double bondi_density(double rinv) { return std::pow(rinv, 1.2); }

/**
 * @brief Get the value of the Bondi velocity at the given inverse radius.
 *
 * @param rinv Inverse radius (in internal units of L^-1).
 * @return Value of the fluid velocity (in internal units of L T^-1).
 */
double bondi_velocity(double rinv) { return -std::pow(rinv, 0.8); }

/**
 * @brief Get the value of the Bondi pressure at the given inverse radius.
 *
 * @param rinv Inverse radius (in internal units of L^-1).
 * @return Value of the pressure (in internal units of M L^-1 T^-2).
 */
double bondi_pressure(double rinv) {
  return ISOTHERMAL_C_SQUARED * bondi_density(rinv);
}

#if IC == IC_BONDI
/**
 * @brief Initialize the given cells.
 *
 * @param cells Cells to initialize.
 */
inline static void initialize(Cell cells[NCELL + 2]) {

  for (unsigned int i = 1; i < NCELL + 1; ++i) {
    const double r_inv = 0.9 / cells[i]._midpoint;
    cells[i]._rho = bondi_density(r_inv);
    cells[i]._u = bondi_velocity(r_inv);
    cells[i]._P = bondi_pressure(r_inv);
    const double r2 = cells[i]._midpoint * cells[i]._midpoint;
    cells[i]._a = -G * MASS_POINT_MASS / r2;
    cells[i]._nfac = 0.;
  }
}
#endif

#endif // BONDI_HPP
