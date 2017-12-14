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
 * @file Cell.hpp
 *
 * @brief Single cell of the grid.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef CELL_HPP
#define CELL_HPP

#include <cstdint>

/**
 * @brief Single cell of the grid.
 */
class Cell {
public:
  // primitive (volume dependent) variables

  /*! @brief Density (in internal units of M L^-3). */
  double _rho;
  /*! @brief Fluid velocity (in internal units of L T^-1). */
  double _u;
  /*! @brief Pressure (in internal units of M L^-1 T^-2). */
  double _P;

  /*! @brief Density gradient (in internal units of M L^-4). */
  double _grad_rho;
  /*! @brief Velocity gradient (in internal units of T^-1). */
  double _grad_u;
  /*! @brief Pressure gradient (in internal units of M L^-2 T^-2). */
  double _grad_P;

  // conserved variables

  /*! @brief Mass (in internal units of M). */
  double _m;
  /*! @brief Momentum (in internal units of L M T^-1). */
  double _p;
  /*! @brief Total energy (in internal units of M L^2 T^-2). */
  double _E;

  // geometrical quantities

  /*! @brief 1D cell volume: radial length of the spherical shell (in internal
   *  units of L). */
  double _V;
  /*! @brief Radial coordinate of the center of the shell (in internal units of
   *  L). */
  double _midpoint;
  /*! @brief Radial coordinate of the lower boundary of the shell (in internal
   *  units of L). */
  double _lowlim;
  /*! @brief Radial coordinate of the upper boundary of the shell (in internal
   *  units of L). */
  double _uplim;

  // gravitational quantities

  /*! @brief Gravitational acceleration (in internal units of L T^-2). */
  double _a;

  // ionization quantities

  /*! @brief Neutral fraction. */
  double _nfac;

  // time step

  /*! @brief Integer time step (in integer time units). */
  uint_fast64_t _integer_dt;

  /*! @brief Physical time step (in internal units of T). */
  double _dt;

  // ADAPTIVE OUTPUT

  /*! @brief Index of the cell (for identification in the log file). */
  uint_least16_t _index;

  /*! @brief Last log file entry offset. */
  uint_least64_t _last_entry;

  /*! @brief Last density value that was written to the log file (in internal
   *  units of M L^-3). */
  double _last_rho;

  /*! @brief Last velocity value that was written to the log file (in internal
   *  units of L T^-1). */
  double _last_u;

  /*! @brief Last pressure value that was written to the log file (in internal
   *  units of M L^-1 T^-2). */
  double _last_P;

  /*! @brief Last neutral fraction value that was written to the log file. */
  double _last_nfac;
};

#endif // CELL_HPP
