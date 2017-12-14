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
 * @file SafeParameters.hpp
 *
 * @brief File that should be included by all files that need parameters.
 * Contains sanity checks on the chosen parameter values. Parameter values
 * should be set in Parameters.hpp.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef SAFEPARAMETERS_HPP
#define SAFEPARAMETERS_HPP

// first include the option names
#include "OptionNames.hpp"

// now include the actual parameters
#include "Parameters.hpp"

// check the parameters with options

/*! @brief Macro that prints the actual value of a macro. */
#define value_of_macro_as_string(macro) #macro
/*! @brief Macro that prints a macro name and its actual value. */
#define value_of_macro(macro) #macro ": " value_of_macro_as_string(macro)

// check boundary conditions
#ifndef BOUNDARIES
#error "No boundary conditions selected!"
#else
#if BOUNDARIES != BOUNDARIES_OPEN && BOUNDARIES != BOUNDARIES_REFLECTIVE &&    \
    BOUNDARIES != BOUNDARIES_BONDI
#pragma message(value_of_macro(BOUNDARIES))
#error "Invalid boundary conditions selected!"
#endif
#endif

// check the equation of state
#ifndef EOS
#error "No equation of state selected!"
#else
#if EOS != EOS_IDEAL && EOS != EOS_ISOTHERMAL && EOS != EOS_BONDI
#pragma message(value_of_macro(EOS))
#error "Invalid equation of state selected!"
#endif
#endif

// check the potential
#ifndef POTENTIAL
#error "No external potential selected!"
#else
#if POTENTIAL != POTENTIAL_NONE && POTENTIAL != POTENTIAL_POINT_MASS
#pragma message(value_of_macro(POTENTIAL))
#error "Invalid potential selected!"
#endif
#endif

// check the initial condition
#ifndef IC
#error "No initial condition selected!"
#else
#if IC != IC_SOD && IC != IC_BONDI && IC != IC_FILE
#pragma message(value_of_macro(IC))
#error "Invalid initial condition selected!"
#endif
#endif

// check ionisation mode
#ifndef IONISATION_MODE
#error "No ionisation mode selected!"
#else
#if IONISATION_MODE != IONISATION_MODE_SELF_CONSISTENT &&                      \
    IONISATION_MODE != IONISATION_MODE_CONSTANT
#pragma message(value_of_macro(IONISATION_MODE))
#error "Invalid ionisation mode selected!"
#endif
#endif

// check ionisation transition
#ifndef IONISATION_TRANSITION
#error "No ionisation transition type selected!"
#else
#if IONISATION_TRANSITION != IONISATION_TRANSITION_JUMP &&                     \
    IONISATION_TRANSITION != IONISATION_TRANSITION_SMOOTH
#pragma message(value_of_macro(IONISATION_TRANSITION))
#error "Invalid ionisation transition type selected!"
#endif
#endif

// check Riemann solver type
#ifndef RIEMANNSOLVER_TYPE
#error "No Riemann solver selected!"
#else
#if RIEMANNSOLVER_TYPE != RIEMANNSOLVER_TYPE_EXACT &&                          \
    RIEMANNSOLVER_TYPE != RIEMANNSOLVER_TYPE_HLLC
#pragma message(value_of_macro(RIEMANNSOLVER_TYPE))
#error "Invalid Riemann solver selected!"
#endif
#endif

// include derived parameters
#include "DerivedParameters.hpp"

#endif // SAFEPARAMETERS_HPP
