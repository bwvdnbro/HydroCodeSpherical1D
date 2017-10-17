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
 * @file Spherical.hpp
 *
 * @brief Spherical source terms that make the 1D method 3D spherical.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef SPHERICAL_HPP
#define SPHERICAL_HPP

/**
 * @brief Add the spherical source terms.
 *
 * See Toro, 2009, chapter 17.
 * We use a second order Runge-Kutta step and apply an operator splitting method
 * to couple the source term to the hydro step.
 */
#define add_spherical_source_term()                                            \
  for (unsigned int i = 1; i < NCELL + 1; ++i) {                               \
    const double r = cells[i]._midpoint;                                       \
    const double Ui[3] = {cells[i]._m / cells[i]._V,                           \
                          cells[i]._p / cells[i]._V,                           \
                          cells[i]._E / cells[i]._V};                          \
    double K1[3];                                                              \
    K1[0] = -DT * Ui[1] / r;                                                   \
    K1[1] = -DT * Ui[1] * Ui[1] / Ui[0] / r;                                   \
    const double p1 = (GAMMA - 1.) * (Ui[2] - 0.5 * Ui[1] * Ui[1] / Ui[0]);    \
    K1[2] = -DT * Ui[1] * (Ui[2] + p1) / Ui[0] / r;                            \
    double U[3] = {Ui[0] + K1[0], Ui[1] + K1[1], Ui[2] + K1[2]};               \
    double K2[3];                                                              \
    K2[0] = -DT * U[1] / r;                                                    \
    K2[1] = -DT * U[1] * U[1] / U[0] / r;                                      \
    const double p2 = (GAMMA - 1.) * (U[2] - 0.5 * U[1] * U[1] / U[0]);        \
    K2[2] = -DT * U[1] * (U[2] + p2) / U[0] / r;                               \
    U[0] = Ui[0] + 0.5 * (K1[0] + K2[0]);                                      \
    U[1] = Ui[1] + 0.5 * (K1[1] + K2[1]);                                      \
    U[2] = Ui[2] + 0.5 * (K1[2] + K2[2]);                                      \
                                                                               \
    cells[i]._m = U[0] * cells[i]._V;                                          \
    cells[i]._p = U[1] * cells[i]._V;                                          \
    cells[i]._E = U[2] * cells[i]._V;                                          \
  }

#endif // SPHERICAL_HPP
