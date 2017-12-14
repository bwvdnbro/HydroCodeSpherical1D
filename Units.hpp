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
 * @file Units.hpp
 *
 * @brief Computation of units not present in Parameters.hpp.
 *
 * The parameters set the units of length and mass, plus a value for the
 * gravitational constant \f$G\f$. If \f$L_{si}\f$, \f$M_{si}\f$ and
 * \f$T_{si}\f$ are the internal units of length, mass and time (expressed in
 * SI units), and \f$G_{si}\f$ and \f$G_{int}\f$ are the values of the
 * gravitational constant in respectively SI units and internal units, then we
 * have
 * \f[
 *   G_{int} \frac{L_{si}^3}{M_{si}T_{si}^2} = G_{si},
 * \f]
 * or
 * \f[
 *   T_{si} = \sqrt{\frac{G_{int}}{G_{si}} \frac{L_{si}^3}{M_{si}}}.
 * \f]
 * This allows us to compute the time unit $T_{si}$. Once we have the units of
 * the three main quantities, all other hydrodynamical units are trivial to
 * compute.
 *
 * Temperatures are always assumed to be in SI units. Other quantities are not
 * used in this problem.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef UNITS_HPP
#define UNITS_HPP

// Internal units in SI units

/*! @brief Time unit (in s). */
#define UNIT_TIME_IN_SI                                                        \
  std::sqrt(                                                                   \
      G_INTERNAL *UNIT_LENGTH_IN_SI *UNIT_LENGTH_IN_SI *UNIT_LENGTH_IN_SI /    \
      (UNIT_MASS_IN_SI * NEWTON_G_IN_SI))

/*! @brief Density unit (in kg m^-3). */
#define UNIT_DENSITY_IN_SI                                                     \
  (UNIT_MASS_IN_SI /                                                           \
   (UNIT_LENGTH_IN_SI * UNIT_LENGTH_IN_SI * UNIT_LENGTH_IN_SI))

/*! @brief Velocity unit (in m s^-1). */
#define UNIT_VELOCITY_IN_SI (UNIT_LENGTH_IN_SI / UNIT_TIME_IN_SI)

/*! @brief Pressure unit (in kg m^-1 s^-2). */
#define UNIT_PRESSURE_IN_SI                                                    \
  (UNIT_MASS_IN_SI / (UNIT_LENGTH_IN_SI * UNIT_TIME_IN_SI * UNIT_TIME_IN_SI))

// Non SI unit conversions

/*! @brief Mass unit (in Msol). */
#define UNIT_MASS_IN_MSOL (UNIT_MASS_IN_SI / SOLAR_MASS_IN_SI)

/*! @brief Time unit (in yr). */
#define UNIT_TIME_IN_YR (UNIT_TIME_IN_SI / YEAR_IN_SI)

/*! @brief Lenght unit (in AU). */
#define UNIT_LENGTH_IN_AU (UNIT_LENGTH_IN_SI / AU_IN_SI)

#endif // UNITS_HPP
