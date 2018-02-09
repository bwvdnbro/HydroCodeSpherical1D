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
 * @file IC.hpp
 *
 * @brief Set up initial conditions.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef IC_HPP
#define IC_HPP

#include "SafeParameters.hpp"

#if IC == IC_BONDI
#include "Bondi.hpp"
#elif IC == IC_SOD
#include "Sod.hpp"
#elif IC == IC_FILE
#include "ICFile.hpp"
#elif IC == IC_BLASTWAVES
#include "BlastWaves.hpp"
#endif

#endif // IC_HPP
