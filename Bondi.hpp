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
#include "SafeParameters.hpp"

#include <cmath>

#if EOS == EOS_BONDI

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
  double Cion =                                                                \
      5. * std::pow(0.9, 2.4) / 3. *                                           \
      (std::pow(INITIAL_IONIZATION_RADIUS, 0.6) - std::pow(0.1, 0.6));         \
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
  for (unsigned int i = 1; i < NCELL + 1; ++i) {                               \
    cells[i]._nfac = get_neutral_fraction(cells[i]._midpoint, rion);           \
  }                                                                            \
  for (unsigned int i = 1; i < NCELL + 1; ++i) {                               \
    const double nfac = cells[i]._nfac;                                        \
    const double ifac = 1. - nfac;                                             \
    cells[i]._P = ISOTHERMAL_C_SQUARED * cells[i]._rho * (100. * ifac + nfac); \
  }

#endif

#if BOUNDARIES == BOUNDARIES_BONDI

/**
 * @brief Apply boundary conditions after the primitive variable conversion.
 */
#define boundary_conditions_primitive_variables()                              \
  /* impose the Bondi solution at the boundaries */                            \
  const double r_inv_low = 0.9 / cells[0]._midpoint;                           \
  cells[0]._rho = bondi_density(r_inv_low);                                    \
  cells[0]._u = bondi_velocity(r_inv_low);                                     \
  if (istep < NSTEP_RELAX) {                                                   \
    cells[0]._P = bondi_pressure(r_inv_low);                                   \
  } else {                                                                     \
    cells[0]._P = 100. * bondi_pressure(r_inv_low);                            \
  }                                                                            \
                                                                               \
  const double r_inv_high = 0.9 / cells[NCELL + 1]._midpoint;                  \
  cells[NCELL + 1]._rho = bondi_density(r_inv_high);                           \
  cells[NCELL + 1]._u = bondi_velocity(r_inv_high);                            \
  cells[NCELL + 1]._P = bondi_pressure(r_inv_high);

/**
 * @brief Apply boundary conditions after the gradient computation.
 */
#define boundary_conditions_gradients()                                        \
  /* compute the exact value of the gradient for the bondi solution */         \
  const double rmin = cells[0]._midpoint - CELLSIZE;                           \
  const double rmax = cells[NCELL + 1]._midpoint + CELLSIZE;                   \
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
      Pmin = cells[0]._P - 100. * bondi_pressure(rmin_inv);                    \
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
    rhomin = cells[NCELL + 1]._rho - cells[NCELL]._rho;                        \
    rhoplu = cells[NCELL + 1]._rho - bondi_density(rmax_inv);                  \
    umin = cells[NCELL + 1]._u - cells[NCELL]._u;                              \
    uplu = cells[NCELL + 1]._u - bondi_velocity(rmax_inv);                     \
    Pmin = cells[NCELL + 1]._P - cells[NCELL]._P;                              \
    Pplu = cells[NCELL + 1]._P - bondi_pressure(rmax_inv);                     \
    if (std::abs(rhomin) < std::abs(rhoplu)) {                                 \
      cells[NCELL + 1]._grad_rho = rhomin / CELLSIZE;                          \
    } else {                                                                   \
      cells[NCELL + 1]._grad_rho = rhoplu / CELLSIZE;                          \
    }                                                                          \
    if (std::abs(umin) < std::abs(uplu)) {                                     \
      cells[NCELL + 1]._grad_u = umin / CELLSIZE;                              \
    } else {                                                                   \
      cells[NCELL + 1]._grad_u = uplu / CELLSIZE;                              \
    }                                                                          \
    if (std::abs(Pmin) < std::abs(Pplu)) {                                     \
      cells[NCELL + 1]._grad_P = Pmin / CELLSIZE;                              \
    } else {                                                                   \
      cells[NCELL + 1]._grad_P = Pplu / CELLSIZE;                              \
    }                                                                          \
  }

#endif

/**
 * @brief Get the value of the Bondi density at the given inverse radius.
 *
 * @param rinv Inverse radius (in internal units of L^-1).
 * @return Value of the density (in internal units of M L^-3).
 */
double bondi_density(double rinv) { return std::pow(rinv, 1.2); }

/**
 * @brief Get the value of the Bondi velocity at the given inverse radius.
 *
 * @param rinv Inverse radius (in internal units of L^-1).
 * @return Value of the fluid velocity (in internal units of L T^-1).
 */
double bondi_velocity(double rinv) { return -std::pow(rinv, 0.8); }

/**
 * @brief Get the value of the Bondi pressure at the given inverse radius.
 *
 * @param rinv Inverse radius (in internal units of L^-1).
 * @return Value of the pressure (in internal units of M L^-1 T^-2).
 */
double bondi_pressure(double rinv) {
  return ISOTHERMAL_C_SQUARED * bondi_density(rinv);
}

#if IC == IC_BONDI
/**
 * @brief Initialize the given cells.
 *
 * @param cells Cells to initialize.
 */
inline static void initialize(Cell cells[NCELL + 2]) {

  for (unsigned int i = 1; i < NCELL + 1; ++i) {
    const double r_inv = 0.9 / cells[i]._midpoint;
    cells[i]._rho = bondi_density(r_inv);
    cells[i]._u = bondi_velocity(r_inv);
    cells[i]._P = bondi_pressure(r_inv);
    const double r2 = cells[i]._midpoint * cells[i]._midpoint;
    cells[i]._a = -G * MASS_POINT_MASS / r2;
    cells[i]._nfac = 0.;
  }
}
#endif

#endif // BONDI_HPP
