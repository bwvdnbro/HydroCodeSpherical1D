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
 * @file EOS.hpp
 *
 * @brief Equation of state.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef EOS_HPP
#define EOS_HPP

// The Bondi specific EOS stuff is in the file Bondi.hpp
#if EOS == EOS_BONDI
#include "Bondi.hpp"
#else

/**
 * @brief Set the initial value for the pressure of the given cell.
 *
 * @param cell Cell.
 */
#if EOS == EOS_IDEAL
#define initial_pressure(cell)
#elif EOS == EOS_ISOTHERMAL
#define initial_pressure(cell) cell._P = ISOTHERMAL_C_SQUARED * cell._rho
#endif

/**
 * @brief Code to determine the neutral fraction of the cells.
 */
#define do_ionization()

/**
 * @brief Conversion function called during the primitive variable conversion
 * for the given cell.
 *
 * @param cell Cell.
 */
#if EOS == EOS_IDEAL
#define update_pressure(cell)                                                  \
  cell._P = (GAMMA - 1.) *                                                     \
            (cell._E / cell._V - 0.5 * cell._rho * cell._u * cell._u);
#elif EOS == EOS_ISOTHERMAL
#define update_pressure(cell) cell._P = ISOTHERMAL_C_SQUARED * cell._rho;
#endif

#endif

#endif // EOS_HPP
