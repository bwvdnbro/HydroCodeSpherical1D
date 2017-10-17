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
 * @file ICFile.hpp
 *
 * @brief Set up initial conditions from a binary initial condition file.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef ICFILE_HPP
#define ICFILE_HPP

#include <fstream>

/**
 * @brief Initialize the given cells.
 *
 * @param cells Cells to initialize.
 */
inline static void initialize(Cell cells[NCELL + 2]) {

  // open the initial condition file
  std::ifstream icfile(IC_FILE_NAME);

  // initialize the cells (we don't initialize the ghost cells)
  for (unsigned int i = 1; i < NCELL + 1; ++i) {
    // read the primitive variables from the initial condition file
    icfile.read(reinterpret_cast<char *>(&cells[i]._rho), sizeof(double));
    icfile.read(reinterpret_cast<char *>(&cells[i]._u), sizeof(double));
    icfile.read(reinterpret_cast<char *>(&cells[i]._P), sizeof(double));
    icfile.read(reinterpret_cast<char *>(&cells[i]._a), sizeof(double));
  }

  // close the initial condition file
  icfile.close();
}

#endif // ICFILE_HPP
