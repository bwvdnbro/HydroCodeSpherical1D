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
 * @file Potential.hpp
 *
 * @brief External potential.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef POTENTIAL_HPP
#define POTENTIAL_HPP

/**
 * @brief Add the gravitational acceleration.
 */
#if POTENTIAL == POTENTIAL_POINT_MASS
#define do_gravity() /* add gravitational acceleration */                      \
  _Pragma("omp parallel for") for (uint_fast32_t i = 1; i < ncell + 1; ++i) {  \
    const double m = cells[i]._V * cells[i]._rho;                              \
    cells[i]._p += 0.5 * cells[i]._dt * cells[i]._a * m;                       \
    /* we do not update the total energy, as we only run gravity simulations   \
       with an isothermal eos, in which case the total energy is ignored by    \
       the hydro scheme */                                                     \
    /*const double r = cells[i]._midpoint;                                     \
    const double a = -G * MASS_POINT_MASS / (r * r);                           \
    cells[i]._a = a;                                                           \
    const double m = cells[i]._V * cells[i]._rho;                              \
    cells[i]._p += 0.5 * DT * a * m; \*/                                       \
  }
#elif POTENTIAL == POTENTIAL_NONE
#define do_gravity()
#endif

/**
 * @brief Do the gravitational half time step prediction for the given cell.
 *
 * @param cell Cell.
 * @param half_dt Half the particle time step (in internal units of T).
 */
#if POTENTIAL != POTENTIAL_NONE
#define add_gravitational_prediction(cell, half_dt)                            \
  cell._u += half_dt * cell._a;
#else
#define add_gravitational_prediction(cell, half_dt)
#endif

#endif // POTENTIAL_HPP
