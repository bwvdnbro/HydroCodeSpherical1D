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
 * @file PhysicalConstants.hpp
 *
 * @brief Values of some physical constants.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef PHYSICALCONSTANTS_HPP
#define PHYSICALCONSTANTS_HPP

/*! @brief Newton gravity constant G (in m^3 kg^-1 s^-2). */
#define NEWTON_G_IN_SI 6.67408e-11

/*! @brief Boltzmann constant k (in m^2 kg s^-1 K^-1). */
#define BOLTZMANN_K_IN_SI 1.38064852e-23

/*! @brief Hydrogen mass (in kg). */
#define HYDROGEN_MASS_IN_SI 1.674e-27

/*! @brief Solar mass (in kg). */
#define SOLAR_MASS_IN_SI 1.9891e30

/*! @brief Year (in s). */
#define YEAR_IN_SI (365.25 * 24. * 3600.)

/*! @brief Astronomical unit (in m). */
#define AU_IN_SI 1.496e11

#endif // PHYSICALCONSTANTS_HPP
