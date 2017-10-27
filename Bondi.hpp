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
 * @file Bondi.hpp
 *
 * @brief Bondi solution expressions.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef BONDI_HPP
#define BONDI_HPP

#include "Cell.hpp"
#include "LambertW.hpp"
#include "SafeParameters.hpp"

#include <cmath>

#if EOS == EOS_BONDI

#define BONDI_DENSITY (BONDI_DENSITY_IN_SI / UNIT_DENSITY_IN_SI)

/*! @brief Bondi radius (in internal units of L). */
#define RBONDI (0.5 * G * MASS_POINT_MASS / ISOTHERMAL_C_SQUARED)

/*! @brief Bondi central density squared (in internal units of M^2 L^-6). */
#define BONDI_DENSITY2 (BONDI_DENSITY * BONDI_DENSITY)

/*! @brief \f$3-2E\f$, with \f$E\f$ the Bondi density exponent. */
#define BONDI_THREE_MINUS_TWO_BONDI_DENSITY_EXPONENT                           \
  (3. - 2. * BONDI_DENSITY_EXPONENT)

/*! @brief Total ionizing luminosity outside the mask, assuming an initial
 *  ionization radius with a given value, and assuming the initial density
 *  profile is the Bondi profile with the given exponent and central density (in
 *  internal units of M^2 L^-3). */
#define BONDI_Q (bondi_Q(INITIAL_IONIZATION_RADIUS) - bondi_Q(RMIN))

/**
 * @brief Set the initial value for the pressure of the given cell.
 *
 * @param cell Cell.
 */
#define initial_pressure(cell) cell._P = ISOTHERMAL_C_SQUARED * cell._rho

/**
 * @brief Code executed before the primitive variable conversion loop is
 * entered.
 */
#define before_primitive_variable_conversion()                                 \
  /* total ionizing budget of the central source                               \
     every shell will absorb a fraction of this budget until the ionization    \
     radius is reached */                                                      \
  double Cion = const_bondi_Q;                                                 \
  double rion = 0.;

/**
 * @brief Conversion function called during the primitive variable conversion
 * for the given cell.
 *
 * @param cell Cell.
 */
#define update_pressure(cell)                                                  \
  /* if ionization is active: check if the cell is ionized */                  \
  if (istep > NSTEP_RELAX && Cion >= 0.) {                                     \
    /* as long as there is ionizing radiation left we assume cells are         \
       ionized */                                                              \
    const double rmin = cells[i]._midpoint - HALF_CELLSIZE;                    \
    const double rmax = cells[i]._midpoint + HALF_CELLSIZE;                    \
    const double Vshell = (rmax * rmax * rmax - rmin * rmin * rmin) / 3.;      \
    const double Cshell = Vshell * cells[i]._rho * cells[i]._rho;              \
    const double ifac = std::min(1., Cion / Cshell);                           \
    if (ifac < 1.) {                                                           \
      if (ifac > 0.) {                                                         \
        const double nfac = 1. - ifac;                                         \
        rion = ifac * rmax + nfac * rmin;                                      \
      }                                                                        \
    } else {                                                                   \
      rion = rmax;                                                             \
    }                                                                          \
    /* subtract this shell's ionization budget from the total */               \
    Cion -= ifac * Cshell;                                                     \
  }

/**
 * @brief Code executed after the primitive variable conversion loop is
 * executed.
 */
#define after_primitive_variable_conversion()                                  \
  for (unsigned int i = 1; i < ncell + 1; ++i) {                               \
    cells[i]._nfac = get_neutral_fraction(cells[i]._midpoint - HALF_CELLSIZE,  \
                                          cells[i]._midpoint + HALF_CELLSIZE,  \
                                          rion, transition_width);             \
  }                                                                            \
  for (unsigned int i = 1; i < ncell + 1; ++i) {                               \
    const double nfac = cells[i]._nfac;                                        \
    const double ifac = 1. - nfac;                                             \
    cells[i]._P = ISOTHERMAL_C_SQUARED * cells[i]._rho *                       \
                  (bondi_pressure_contrast * ifac + nfac);                     \
  }

#endif

#if BOUNDARIES == BOUNDARIES_BONDI

/**
 * @brief Apply boundary conditions after the primitive variable conversion.
 */
#define boundary_conditions_primitive_variables()                              \
  /* impose the Bondi solution at the boundaries */                            \
  const double r_inv_low = RBONDI / cells[0]._midpoint;                        \
  cells[0]._rho = bondi_density(r_inv_low);                                    \
  cells[0]._u = bondi_velocity(r_inv_low);                                     \
  if (istep < NSTEP_RELAX) {                                                   \
    cells[0]._P = bondi_pressure(r_inv_low);                                   \
  } else {                                                                     \
    cells[0]._P = bondi_pressure_contrast * bondi_pressure(r_inv_low);         \
  }                                                                            \
                                                                               \
  const double r_inv_high = RBONDI / cells[ncell + 1]._midpoint;               \
  cells[ncell + 1]._rho = bondi_density(r_inv_high);                           \
  cells[ncell + 1]._u = bondi_velocity(r_inv_high);                            \
  cells[ncell + 1]._P = bondi_pressure(r_inv_high);

/**
 * @brief Apply boundary conditions after the gradient computation.
 */
#define boundary_conditions_gradients()                                        \
  /* compute the exact value of the gradient for the bondi solution */         \
  const double rmin = cells[0]._midpoint - CELLSIZE;                           \
  const double rmax = cells[ncell + 1]._midpoint + CELLSIZE;                   \
  const double rmin_inv = 1. / rmin;                                           \
  const double rmax_inv = 1. / rmax;                                           \
                                                                               \
  /* lower boundary */                                                         \
  {                                                                            \
    double rhomin, rhoplu, umin, uplu, Pmin, Pplu, dmin, dplu;                 \
    rhomin = cells[0]._rho - bondi_density(rmin_inv);                          \
    rhoplu = cells[0]._rho - cells[1]._rho;                                    \
    umin = cells[0]._u - bondi_velocity(rmin_inv);                             \
    uplu = cells[0]._u - cells[1]._u;                                          \
    if (istep < NSTEP_RELAX) {                                                 \
      Pmin = cells[0]._P - bondi_pressure(rmin_inv);                           \
    } else {                                                                   \
      Pmin = cells[0]._P - bondi_pressure_contrast * bondi_pressure(rmin_inv); \
    }                                                                          \
    Pplu = cells[0]._P - cells[1]._P;                                          \
    if (std::abs(rhomin) < std::abs(rhoplu)) {                                 \
      cells[0]._grad_rho = rhomin / CELLSIZE;                                  \
    } else {                                                                   \
      cells[0]._grad_rho = rhoplu / CELLSIZE;                                  \
    }                                                                          \
    if (std::abs(umin) < std::abs(uplu)) {                                     \
      cells[0]._grad_u = umin / CELLSIZE;                                      \
    } else {                                                                   \
      cells[0]._grad_u = uplu / CELLSIZE;                                      \
    }                                                                          \
    if (std::abs(Pmin) < std::abs(Pplu)) {                                     \
      cells[0]._grad_P = Pmin / CELLSIZE;                                      \
    } else {                                                                   \
      cells[0]._grad_P = Pplu / CELLSIZE;                                      \
    }                                                                          \
  }                                                                            \
                                                                               \
  /* upper boundary */                                                         \
  {                                                                            \
    double rhomin, rhoplu, umin, uplu, Pmin, Pplu, dmin, dplu;                 \
    rhomin = cells[ncell + 1]._rho - cells[ncell]._rho;                        \
    rhoplu = cells[ncell + 1]._rho - bondi_density(rmax_inv);                  \
    umin = cells[ncell + 1]._u - cells[ncell]._u;                              \
    uplu = cells[ncell + 1]._u - bondi_velocity(rmax_inv);                     \
    Pmin = cells[ncell + 1]._P - cells[ncell]._P;                              \
    Pplu = cells[ncell + 1]._P - bondi_pressure(rmax_inv);                     \
    if (std::abs(rhomin) < std::abs(rhoplu)) {                                 \
      cells[ncell + 1]._grad_rho = rhomin / CELLSIZE;                          \
    } else {                                                                   \
      cells[ncell + 1]._grad_rho = rhoplu / CELLSIZE;                          \
    }                                                                          \
    if (std::abs(umin) < std::abs(uplu)) {                                     \
      cells[ncell + 1]._grad_u = umin / CELLSIZE;                              \
    } else {                                                                   \
      cells[ncell + 1]._grad_u = uplu / CELLSIZE;                              \
    }                                                                          \
    if (std::abs(Pmin) < std::abs(Pplu)) {                                     \
      cells[ncell + 1]._grad_P = Pmin / CELLSIZE;                              \
    } else {                                                                   \
      cells[ncell + 1]._grad_P = Pplu / CELLSIZE;                              \
    }                                                                          \
  }

#endif

/**
 * @brief Squared Bondi velocity divided by the sound speed squared.
 *
 * @param rinv Inverse radius (in units of RBONDI^-1).
 * @return Bondi velocity squared divided by the sound speed squared.
 */
double u2_over_cs2(double rinv) {
  const double lambertarg = -std::exp(3. + 4. * (std::log(rinv) - rinv));
  if (rinv < 1.) {
    return -LambertW::lambert_w(lambertarg, 0);
  } else {
    return -LambertW::lambert_w(lambertarg, -1);
  }
}

/**
 * @brief Get the value of the Bondi density at the given inverse radius.
 *
 * @param rinv Inverse radius (in units of RBONDI^-1).
 * @return Value of the density (in internal units of M L^-3).
 */
double bondi_density(double rinv) {
  // we need to manually disable the density very close to r = 0 to prevent
  // errors in the Lambert W function
  if (rinv < 150.) {
    return BONDI_DENSITY * std::exp(-0.5 * u2_over_cs2(rinv) + 2. * rinv - 1.5);
  } else {
    return 0.;
  }
}

/**
 * @brief Get the value of the Bondi velocity at the given inverse radius.
 *
 * @param rinv Inverse radius (in internal units of L^-1).
 * @return Value of the fluid velocity (in internal units of L T^-1).
 */
double bondi_velocity(double rinv) {
  return -std::sqrt(ISOTHERMAL_C_SQUARED * u2_over_cs2(rinv));
}

/**
 * @brief Get the value of the Bondi pressure at the given inverse radius.
 *
 * @param rinv Inverse radius (in internal units of L^-1).
 * @return Value of the pressure (in internal units of M L^-1 T^-2).
 */
double bondi_pressure(double rinv) {
  return ISOTHERMAL_C_SQUARED * bondi_density(rinv);
}

/**
 * @brief Integrand for the luminosity integral.
 *
 * @param r Radius (in internal units of L).
 * @return Integrand: \f$r^2 \rho{}(r)^2\f$ (in internal units of M^2 L^-4).
 */
double bondi_Q_integrand(double r) {
  const double rinv = RBONDI / r;
  const double rho = bondi_density(rinv);
  return r * r * rho * rho;
}

/**
 * @brief Do a simple first order integration of the bondi_Q_integrand over the
 * given interval.
 *
 * @param a Lower limit of the interval (in internal units of L).
 * @param b Upper limit of the interval (in internal units of L).
 * @return First order approximation to the integral over the interval.
 */
double bondi_Q_interval(double a, double b) {
  return (b - a) * bondi_Q_integrand(0.5 * (a + b));
}

/**
 * @brief Get the total luminosity absorbed by the volume within the given
 * radius for the bondi density profile.
 *
 * @param r Radius (in internal units of L).
 * @param tolerance Required relative accuracy for the result (default: 1.e-10).
 * @return Total luminosity absorbed by the volume within the given radius
 * (in internal units of M^2 L^-3).
 */
double bondi_Q(double r, double tolerance = 1.e-8) {
  double I0 = bondi_Q_interval(0., r);
  double I1 = bondi_Q_interval(0., 0.5 * r) + bondi_Q_interval(0.5 * r, r);
  unsigned int npiece = 2;
  while (std::abs(I0 - I1) > std::abs(I0 + I1) * tolerance) {
    npiece <<= 1;
    I0 = I1;
    I1 = 0.;
    const double dr = r / npiece;
    for (unsigned int i = 0; i < npiece; ++i) {
      I1 += bondi_Q_interval(i * dr, (i + 1) * dr);
    }
  }
  return I1;
}

#if IC == IC_BONDI
/**
 * @brief Initialize the given cells.
 *
 * @param cells Cells to initialize.
 */
#define initialize(cells, ncell)                                               \
  for (unsigned int i = 1; i < ncell + 1; ++i) {                               \
    const double r_inv = RBONDI / cells[i]._midpoint;                          \
    cells[i]._rho = bondi_density(r_inv);                                      \
    cells[i]._u = bondi_velocity(r_inv);                                       \
    cells[i]._P = bondi_pressure(r_inv);                                       \
    const double r2 = cells[i]._midpoint * cells[i]._midpoint;                 \
    cells[i]._a = -G * MASS_POINT_MASS / r2;                                   \
    cells[i]._nfac = 0.;                                                       \
  }

#endif

#endif // BONDI_HPP
