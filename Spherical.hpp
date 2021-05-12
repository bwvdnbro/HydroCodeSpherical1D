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
 * @file Spherical.hpp
 *
 * @brief Spherical source terms that make the 1D method 3D spherical.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef SPHERICAL_HPP
#define SPHERICAL_HPP

/**
 * @brief Compute the spherical source term for the given specific conserved
 * variables U (conserved quantities divided by cell volume).
 *
 * @param fac Pre-factor to multiply in (alpha/r * dt).
 * @param U Specific conserved variables.
 * @param S Output source terms.
 */
static inline void compute_source_term(const double fac, const double U[3],
                                       double S[3]) {
  const double U0inv = 1. / U[0];
  const double U12 = U[1] * U[1];
  const double p = (GAMMA - 1.) * (U[2] - 0.5 * U12 * U0inv);
  S[0] = fac * U[1];
  S[1] = fac * U12 * U0inv;
  S[2] = fac * U[1] * U0inv * (U[2] + p);
}

/**
 * @brief Check if the explicit time step satisfies the stability criterion.
 *
 * @param dt Time step.
 * @param Sfac Pre-factor for the source terms (alpha / r).
 * @param Specific conserved variables.
 */
static inline void check_explicit_timestep(const double dt, const double Sfac,
                                           const double U[3]) {
  if (U[1] != 0.) {
    const double maxdt = U[0] / (GAMMA * Sfac * std::abs(U[1]));
    if (dt >= maxdt) {
      std::cerr << "Unstable solution for explicit geometric source term!"
                << std::endl;
      exit(1);
    }
  }
}

/**
 * @brief Implicitly compute new specific conserved variables based on the given
 * initial values, time step and source term pre-factor.
 *
 * This version uses a second order Runge-Kutta scheme.
 *
 * Eq (15.36) in Toro (2009).
 *
 * @param dt Time step.
 * @param Sfac Pre-factor for the source terms (alpha / r).
 * @param Ui Initial values of the specific conserved variables.
 * @param U New values for the specific conserved variables (output).
 */
static inline void second_order_runge_kutta(const double dt, const double Sfac,
                                            const double Ui[3], double U[3]) {

  check_explicit_timestep(dt, Sfac, Ui);

  double K1[3], K2[3];
  compute_source_term(-dt * Sfac, Ui, K1);
  U[0] = Ui[0] + K1[0];
  U[1] = Ui[1] + K1[1];
  U[2] = Ui[2] + K1[2];
  compute_source_term(-dt * Sfac, U, K2);
  U[0] = Ui[0] + 0.5 * (K1[0] + K2[0]);
  U[1] = Ui[1] + 0.5 * (K1[1] + K2[1]);
  U[2] = Ui[2] + 0.5 * (K1[2] + K2[2]);
}

/**
 * @brief Implicitly compute new specific conserved variables based on the given
 * initial values, time step and source term pre-factor.
 *
 * This version uses a fourth order Runge-Kutta scheme.
 *
 * Eq (15.37) in Toro (2009).
 *
 * @param dt Time step.
 * @param Sfac Pre-factor for the source terms (alpha / r).
 * @param Ui Initial values of the specific conserved variables.
 * @param U New values for the specific conserved variables (output).
 */
static inline void fourth_order_runge_kutta(const double dt, const double Sfac,
                                            const double Ui[3], double U[3]) {

  check_explicit_timestep(dt, Sfac, Ui);

  double K1[3], K2[3], K3[3], K4[3];
  compute_source_term(-dt * Sfac, Ui, K1);
  U[0] = Ui[0] + 0.5 * K1[0];
  U[1] = Ui[1] + 0.5 * K1[1];
  U[2] = Ui[2] + 0.5 * K1[2];
  compute_source_term(-dt * Sfac, U, K2);
  U[0] = Ui[0] + 0.5 * K2[0];
  U[1] = Ui[1] + 0.5 * K2[1];
  U[2] = Ui[2] + 0.5 * K2[2];
  compute_source_term(-dt * Sfac, U, K3);
  U[0] = Ui[0] + K3[0];
  U[1] = Ui[1] + K3[1];
  U[2] = Ui[2] + K3[2];
  compute_source_term(-dt * Sfac, U, K4);
  U[0] = Ui[0] + (K1[0] + 2. * K2[0] + 2. * K3[0] + K4[0]) / 6.;
  U[1] = Ui[1] + (K1[1] + 2. * K2[1] + 2. * K3[1] + K4[1]) / 6.;
  U[2] = Ui[2] + (K1[2] + 2. * K2[2] + 2. * K3[2] + K4[2]) / 6.;
}

/**
 * @brief Solve implicitly for the geometrical source terms using the second
 * order accurate trapezoidal method.
 *
 * The implicit method is given by equation (15.35) in Toro (2009):
 * @f[
 *   U^{n+1} = U^n + \frac{1}{2} \Delta{}t [S(U^n) + S(U^{n+1})],
 * @f]
 * with @f$S(U)@f$ defined in compute_source_term().
 *
 * By explicitly writing out the expressions for @f$U=\rho{}, \rho{}u, e@f$,
 * it is clear that @f$u_1^{n+1}@f$ can be implicitly computed from the initial
 * state, since it does not depend on @f$u_0^{n+1}@f$ and @f$u_2^{n+1}@f$. Once
 * @f$u_1^{n+1}@f$ is known, it can be used to compute @f$u_0^{n+1}@f$, which
 * again does not depend on @f$u_2^{n+1}@f$. Once we have computed 2 out of 3
 * unknowns, computing the third one is then trivial (yet involves the most
 * complicated expression).
 *
 * @param dt Time step.
 * @param Sfac Pre-factor for the source terms (alpha / r).
 * @param Ui Initial values of the specific conserved variables.
 * @param U New values for the specific conserved variables (output).
 */
static inline void implicit_trapezoidal_method(const double dt,
                                               const double Sfac,
                                               const double Ui[3],
                                               double U[3]) {

  const double A = -0.5 * dt * Sfac;

  // we first solve for the new momentum
  U[1] = Ui[1] * (Ui[0] + A * Ui[1]) / (Ui[0] - A * Ui[1]);
  // now we can solve for the new density
  U[0] = Ui[0] + A * (Ui[1] + U[1]);
  // and finally we can solve for the new energy
  const double U02 = U[0] * U[0];
  const double Ui02 = Ui[0] * Ui[0];
  U[2] = (Ui[2] * Ui[0] * U02 * (Ui[0] + GAMMA * A * Ui[1]) -
          0.5 * (GAMMA - 1.) * A *
              (U02 * Ui[1] * Ui[1] * Ui[1] + Ui02 * U[1] * U[1] * U[1])) /
         (Ui02 * U[0] * (U[0] - GAMMA * A * U[1]));
}

/**
 * @brief Add the spherical source terms.
 *
 * See Toro, 2009, chapter 17.
 * We use a second order Runge-Kutta step and apply an operator splitting method
 * to couple the source term to the hydro step.
 */
#if DIMENSIONALITY == DIMENSIONALITY_1D
#define add_spherical_source_term()
#elif DIMENSIONALITY == DIMENSIONALITY_3D
#define add_spherical_source_term()                                            \
  _Pragma("omp parallel for") for (uint_fast32_t i = 1; i < ncell + 1; ++i) {  \
    if (cells[i]._m > 0.) {                                                    \
      const double r = cells[i]._midpoint;                                     \
      const double Sfac = 1. / r;                                              \
      const double Vinv = 1. / cells[i]._V;                                    \
      const double dt = cells[i]._dt;                                          \
      const double Ui[3] = {cells[i]._m * Vinv, cells[i]._p * Vinv,            \
                            cells[i]._E * Vinv};                               \
      double U[3];                                                             \
      implicit_trapezoidal_method(dt, Sfac, Ui, U);                            \
      cells[i]._m = U[0] * cells[i]._V;                                        \
      cells[i]._p = U[1] * cells[i]._V;                                        \
      cells[i]._E = U[2] * cells[i]._V;                                        \
    }                                                                          \
  }
#endif

#endif // SPHERICAL_HPP
