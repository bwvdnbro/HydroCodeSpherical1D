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
#include <iostream>

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

  // ionisation quantities

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

  /**
   * @brief Print the contents of the cell for debugging purposes.
   */
  inline void print() {
    std::cerr << "Cell values:" << std::endl;
    std::cerr << "_rho: " << _rho << std::endl;
    std::cerr << "_u: " << _u << std::endl;
    std::cerr << "_P: " << _P << std::endl;
    std::cerr << "_grad_rho: " << _grad_rho << std::endl;
    std::cerr << "_grad_u: " << _grad_u << std::endl;
    std::cerr << "_grad_P: " << _grad_P << std::endl;
    std::cerr << "_m: " << _m << std::endl;
    std::cerr << "_p: " << _p << std::endl;
    std::cerr << "_E: " << _E << std::endl;
    std::cerr << "_V: " << _V << std::endl;
    std::cerr << "_midpoint: " << _midpoint << std::endl;
    std::cerr << "_lowlim: " << _lowlim << std::endl;
    std::cerr << "_uplim: " << _uplim << std::endl;
    std::cerr << "_a: " << _a << std::endl;
    std::cerr << "_nfac: " << _nfac << std::endl;
    std::cerr << "_integer_dt: " << _integer_dt << std::endl;
    std::cerr << "_dt: " << _dt << std::endl;
    std::cerr << "_index: " << _index << std::endl;
    std::cerr << "_last_entry: " << _last_entry << std::endl;
    std::cerr << "_last_rho: " << _last_rho << std::endl;
    std::cerr << "_last_u: " << _last_u << std::endl;
    std::cerr << "_last_P: " << _last_P << std::endl;
    std::cerr << "_last_nfac: " << _last_nfac << std::endl;
  };
};

#endif // CELL_HPP
