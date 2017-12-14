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

#define BONDI_DENSITY (BONDI_DENSITY_IN_SI / UNIT_DENSITY_IN_SI)

/*! @brief Bondi radius (in internal units of L). */
#define RBONDI (0.5 * G * MASS_POINT_MASS / ISOTHERMAL_C_SQUARED)
#define RBONDI_ION                                                             \
  (0.5 * G * MASS_POINT_MASS / (bondi_pressure_contrast * ISOTHERMAL_C_SQUARED))

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

#if EOS == EOS_BONDI

/**
 * @brief Initialize ionization variables.
 */
#define ionization_initialize()                                                \
  double rion_old = 0.;                                                        \
                                                                               \
  std::ofstream bondi_rfile("ionization_radius.dat");                          \
  /*bondi_rfile << "# time (s)\tionization radius (m)\n";*/                    \
                                                                               \
  const double bondi_S =                                                       \
      (transition_width > 0.) ? 3. / (4. * transition_width) : 0.;             \
  const double bondi_A = -32. * bondi_S * bondi_S * bondi_S / 27.;             \
                                                                               \
  const double bondi_rmin = RMIN - CELLSIZE;                                   \
  /*const double bondi_volume_correction_factor = 4. * M_PI / 3. / CELLSIZE *  \
    (RMIN * RMIN * RMIN - bondi_rmin * bondi_rmin * bondi_rmin);*/             \
  const double bondi_volume_correction_factor = 0.;                            \
  std::cout << "Bondi volume correction factor: "                              \
            << bondi_volume_correction_factor << std::endl;                    \
                                                                               \
  /*const double const_bondi_Q = 1.47132e-21;*/                                \
  std::cout << "Precomputing Bondi luminosity..." << std::endl;                \
                                                                               \
  _Pragma("omp parallel for") for (uint_fast32_t i = 1; i < ncell + 1; ++i) {  \
    const double rmin = cells[i]._lowlim;                                      \
    const double rmax = cells[i]._uplim;                                       \
    const double Vshell = (rmax * rmax * rmax - rmin * rmin * rmin) / 3.;      \
    const double Cshell = Vshell * cells[i]._rho * cells[i]._rho;              \
    cells[i]._nfac = Cshell;                                                   \
  }                                                                            \
  double const_bondi_Q = 0.;                                                   \
  for (uint_fast32_t i = 1; i < ncell + 1; ++i) {                              \
    const double rmin = cells[i]._lowlim;                                      \
    const double rmax = cells[i]._uplim;                                       \
    if (rmax < INITIAL_IONIZATION_RADIUS) {                                    \
      const_bondi_Q += cells[i]._nfac;                                         \
    } else if (rmin < INITIAL_IONIZATION_RADIUS) {                             \
      const double ifac =                                                      \
          (INITIAL_IONIZATION_RADIUS * INITIAL_IONIZATION_RADIUS *             \
               INITIAL_IONIZATION_RADIUS -                                     \
           rmin * rmin * rmin) /                                               \
          (rmax * rmax * rmax - rmin * rmin * rmin);                           \
      const double Cshell = cells[i]._nfac * ifac;                             \
      const_bondi_Q += Cshell;                                                 \
    }                                                                          \
  }                                                                            \
  std::cout << "Bondi Q: " << const_bondi_Q << std::endl;                      \
  double central_mass = MASS_POINT_MASS;

/**
 * @brief Set the initial value for the pressure of the given cell.
 *
 * @param cell Cell.
 */
#define initial_pressure(cell) cell._P = ISOTHERMAL_C_SQUARED * cell._rho

/**
 * @brief Conversion function called during the primitive variable conversion
 * for the given cell.
 *
 * @param cell Cell.
 */
#define update_pressure(cell)                                                  \
  const double nfac = cells[i]._nfac;                                          \
  const double ifac = 1. - nfac;                                               \
  cells[i]._P = ISOTHERMAL_C_SQUARED * cells[i]._rho *                         \
                (bondi_pressure_contrast * ifac + nfac);

/**
 * @brief Code to determine the neutral fraction of the cells.
 */
#define do_ionization()                                                        \
  _Pragma("omp parallel for") for (uint_fast32_t i = 1; i < ncell + 1; ++i) {  \
    const double rmin = cells[i]._lowlim;                                      \
    const double rmax = cells[i]._uplim;                                       \
    const double Vshell = (rmax * rmax * rmax - rmin * rmin * rmin) / 3.;      \
    const double Cshell = Vshell * cells[i]._rho * cells[i]._rho;              \
    cells[i]._nfac = Cshell;                                                   \
  }                                                                            \
                                                                               \
  double Cion =                                                                \
      const_bondi_Q * get_bondi_Q_factor(central_mass / MASS_POINT_MASS);      \
  double rion = 0.;                                                            \
  for (uint_fast32_t i = 1; i < ncell + 1; ++i) {                              \
    if (Cion > 0.) {                                                           \
      const double rmin = cells[i]._lowlim;                                    \
      const double rmax = cells[i]._uplim;                                     \
      const double ifac = std::min(1., Cion / cells[i]._nfac);                 \
      if (ifac < 1.) {                                                         \
        if (ifac > 0.) {                                                       \
          const double nfac = 1. - ifac;                                       \
          rion = std::cbrt(ifac * rmax * rmax * rmax +                         \
                           nfac * rmin * rmin * rmin);                         \
        }                                                                      \
      } else {                                                                 \
        rion = rmax;                                                           \
      }                                                                        \
      /* subtract this shell's ionization budget from the total */             \
      Cion -= ifac * cells[i]._nfac;                                           \
    }                                                                          \
  }                                                                            \
  if (std::abs(rion - rion_old) > 1.e-2 * std::abs(rion + rion_old)) {         \
    const double curtime =                                                     \
        current_integer_time * time_conversion_factor * UNIT_TIME_IN_SI;       \
    const double ionrad = rion * UNIT_LENGTH_IN_SI;                            \
    Cion = const_bondi_Q * get_bondi_Q_factor(central_mass / MASS_POINT_MASS); \
    bondi_rfile.write(reinterpret_cast<const char *>(&curtime),                \
                      sizeof(double));                                         \
    bondi_rfile.write(reinterpret_cast<const char *>(&ionrad),                 \
                      sizeof(double));                                         \
    bondi_rfile.write(reinterpret_cast<const char *>(&Cion), sizeof(double));  \
    bondi_rfile.flush();                                                       \
    rion_old = rion;                                                           \
  }                                                                            \
  /*rion = INITIAL_IONIZATION_RADIUS;*/                                        \
  const double rion_min = rion - 0.5 * transition_width;                       \
  const double rion_max = rion + 0.5 * transition_width;                       \
  _Pragma("omp parallel for") for (uint_fast32_t i = 1; i < ncell + 1; ++i) {  \
    cells[i]._nfac =                                                           \
        get_neutral_fraction(cells[i]._lowlim, cells[i]._uplim, rion,          \
                             rion_min, rion_max, bondi_S, bondi_A);            \
  }

/**
 * @brief Code to handle the mass flux into the inner mask.
 *
 * @param mflux Mass flux into the inner mask.
 */
#define flux_into_inner_mask(mflux)                                            \
  central_mass -= mflux * bondi_volume_correction_factor;
#endif

#if BOUNDARIES == BOUNDARIES_BONDI

/**
 * @brief Initialize variables used for the boundary conditions.
 */
#define boundary_conditions_initialize()                                       \
  const double bondi_r_inv_low = RBONDI_ION / cells[0]._midpoint;              \
  const double bondi_density_low = bondi_density(bondi_r_inv_low);             \
  const double bondi_velocity_low = bondi_velocity(bondi_r_inv_low);           \
  const double bondi_pressure_low = bondi_pressure(bondi_r_inv_low);           \
                                                                               \
  const double bondi_r_inv_high = RBONDI / cells[ncell + 1]._midpoint;         \
  const double bondi_density_high = bondi_density(bondi_r_inv_high);           \
  const double bondi_velocity_high = bondi_velocity(bondi_r_inv_high);         \
  const double bondi_pressure_high = bondi_pressure(bondi_r_inv_high);         \
                                                                               \
  const double bondi_rmin_inv = RBONDI_ION / (cells[0]._midpoint - CELLSIZE);  \
  const double bondi_rmax_inv =                                                \
      RBONDI / (cells[ncell + 1]._midpoint + CELLSIZE);                        \
  const double bondi_density_min = bondi_density(bondi_rmin_inv);              \
  const double bondi_velocity_min = bondi_velocity(bondi_rmin_inv);            \
  const double bondi_pressure_min = bondi_pressure(bondi_rmin_inv);            \
  const double bondi_density_max = bondi_density(bondi_rmax_inv);              \
  const double bondi_velocity_max = bondi_velocity(bondi_rmax_inv);            \
  const double bondi_pressure_max = bondi_pressure(bondi_rmax_inv);

/**
 * @brief Apply boundary conditions after the primitive variable conversion.
 */
#define boundary_conditions_primitive_variables()                              \
  /* impose the Bondi solution at the boundaries */                            \
  /*cells[0]._rho = bondi_density_low;                                         \
  cells[0]._u = bondi_velocity_low;                                            \
  cells[0]._P = bondi_pressure_contrast * bondi_pressure_low;*/                \
  cells[0]._rho = cells[1]._rho;                                               \
  cells[0]._u = cells[1]._u;                                                   \
  cells[0]._P = cells[1]._P;                                                   \
                                                                               \
  /*if(current_integer_time == 0){                                             \
  cells[ncell + 1]._rho = bondi_density_high*(1. + 1.e-5);                     \
  } else {                                                                     \
  cells[ncell + 1]._rho = bondi_density_high;                                  \
  }\*/                                                                         \
  cells[ncell + 1]._rho = bondi_density_high;                                  \
  cells[ncell + 1]._u = bondi_velocity_high;                                   \
  cells[ncell + 1]._P = bondi_pressure_high;

/**
 * @brief Apply boundary conditions after the gradient computation.
 */
#define boundary_conditions_gradients()                                        \
  /* lower boundary */                                                         \
  /*{                                                                          \
    double rhomin, rhoplu, umin, uplu, Pmin, Pplu, dmin, dplu;                 \
    rhomin = cells[0]._rho - bondi_density_min;                                \
    rhoplu = cells[0]._rho - cells[1]._rho;                                    \
    umin = cells[0]._u - bondi_velocity_min;                                   \
    uplu = cells[0]._u - cells[1]._u;                                          \
    Pmin = cells[0]._P - bondi_pressure_contrast * bondi_pressure_min;         \
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
  }*/                                                                          \
  {                                                                            \
    cells[0]._grad_rho = cells[1]._grad_rho;                                   \
    cells[0]._grad_u = cells[1]._grad_u;                                       \
    cells[0]._grad_P = cells[1]._grad_P;                                       \
  }                                                                            \
                                                                               \
  /* upper boundary */                                                         \
  {                                                                            \
    double rhomin, rhoplu, umin, uplu, Pmin, Pplu, dmin, dplu;                 \
    rhomin = cells[ncell + 1]._rho - cells[ncell]._rho;                        \
    rhoplu = cells[ncell + 1]._rho - bondi_density_max;                        \
    umin = cells[ncell + 1]._u - cells[ncell]._u;                              \
    uplu = cells[ncell + 1]._u - bondi_velocity_max;                           \
    Pmin = cells[ncell + 1]._P - cells[ncell]._P;                              \
    Pplu = cells[ncell + 1]._P - bondi_pressure_max;                           \
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

/**
 * @brief Get the increase in luminosity corresponding to the given increase in
 * central mass.
 *
 * @param M Central mass (in units of the initial central mass).
 */
inline static double get_bondi_Q_factor(const double M) {
  return 7.96185873 * (std::pow(M, 2.47692987) - 1.) + 1.;
}

inline static double
get_neutral_fraction_integral(const double A, const double S, const double rion,
                              const double rmin, const double rmax) {
  const double rmaxrel = rmax - rion;
  const double rminrel = rmin - rion;
  const double rdiff = rmax - rmin;
  const double rmaxrel2 = rmaxrel * rmaxrel;
  const double rmaxrel4 = rmaxrel2 * rmaxrel2;
  const double rminrel2 = rminrel * rminrel;
  const double rminrel4 = rminrel2 * rminrel2;
  return A * (rmaxrel4 - rminrel4) + S * (rmaxrel2 - rminrel2) + 0.5 * rdiff;
}

/**
 * @brief Get the neutral fraction for the cell with the given midpoint radius.
 *
 * New version that contains a linear transition from ionized to neutral.
 *
 * @param rmin Radius of the lower wall of the cell (in internal units of L).
 * @param rmax Radius of the upper wall of the cell (in internal units of L).
 * @param rion Ionization radius (in internal units of L).
 * @param transition_width Width of the transition region (in internal units of
 * L).
 * @return Neutral fraction within the cell.
 */
inline static double get_neutral_fraction(const double rmin, const double rmax,
                                          const double rion,
                                          const double rion_min,
                                          const double rion_max, const double S,
                                          const double A) {
  // Note: we assume rmax - rmin << rion_max - rion_min
  if (rmax < rion_min) {
    return 0.;
  } else if (rmin < rion_min && rmax >= rion_min) {
    return get_neutral_fraction_integral(A, S, rion, rion_min, rmax) /
           (rmax - rmin);
  } else if (rmin >= rion_min && rmax <= rion_max) {
    return get_neutral_fraction_integral(A, S, rion, rmin, rmax) /
           (rmax - rmin);
  } else if (rmin < rion_max && rmax > rion_max) {
    return (get_neutral_fraction_integral(A, S, rion, rmin, rion_max) +
            (rmax - rion_max)) /
           (rmax - rmin);
  } else {
    return 1.;
  }
  //    if(rmax < rion){
  //      return 0;
  //    } else if(rmin < rion){
  //      return (rmax - rion) / (rmax - rmin);
  //    } else {
  //      return 1.;
  //    }
}

#if IC == IC_BONDI
/**
 * @brief Initialize the given cells.
 *
 * @param cells Cells to initialize.
 */
#define initialize(cells, ncell)                                               \
  _Pragma("omp parallel for") for (unsigned int i = 1; i < ncell + 1; ++i) {   \
    /*const double r_inv = RBONDI / cells[i]._midpoint;                        \
    cells[i]._rho = bondi_density(r_inv);                                      \
    cells[i]._u = bondi_velocity(r_inv);                                       \
    cells[i]._P = bondi_pressure(r_inv);*/                                     \
    const double r_inv = RBONDI / RMAX;                                        \
    cells[i]._rho = bondi_density(r_inv);                                      \
    cells[i]._u = bondi_velocity(r_inv);                                       \
    cells[i]._P = ISOTHERMAL_C_SQUARED * BONDI_DENSITY;                        \
    const double r2 = cells[i]._midpoint * cells[i]._midpoint;                 \
    cells[i]._a = -G * MASS_POINT_MASS / r2;                                   \
    cells[i]._nfac = 0.;                                                       \
  }

#endif

#endif // BONDI_HPP
