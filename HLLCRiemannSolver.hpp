/*******************************************************************************
 * This file is part of HydroCodeSpherical1D
 * Copyright (C) 2016, 2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * along with Shadowfax. If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/

/**
 * @file HLLCRiemannSolver.hpp
 *
 * @brief HLLC Riemann solver.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be,
 * bv7@st-andrews.ac.uk)
 */
#ifndef HLLCRIEMANNSOLVER_HPP
#define HLLCRIEMANNSOLVER_HPP

#include <algorithm>
#include <cmath>

/**
 * @brief HLLC Riemann solver.
 *
 * Based on Chapter 10 and more specifically the summary in 10.6 of Toro, 2009,
 * Riemann Solvers and Numerical Methods for Fluid Dynamics.
 */
class HLLCRiemannSolver {
private:
  /*! @brief Adiabatic index \f$\gamma{}\f$. */
  double _gamma;

  /*! @brief \f$\frac{2(\gamma{}+1)}{\gamma{}}\f$. */
  double _hgp1dg;

  /*! @brief \f$\frac{1}{\gamma{}-1}\f$. */
  double _odgm1;

public:
  /**
   * @brief Constructor.
   *
   * @param gamma Adiabatic index \f$\gamma{}\f$.
   */
  HLLCRiemannSolver(double gamma = 5. / 3.)
      : _gamma(gamma), _hgp1dg(0.5 * (_gamma + 1.) / _gamma),
        _odgm1(1. / (_gamma - 1.)) {}

  /**
   * @brief Solve the Riemann problem with the given left and right state
   * directly for the flux.
   *
   * @param rhoL Left state density.
   * @param uL Left state velocity.
   * @param PL Left state pressure.
   * @param rhoR Right state density.
   * @param uR Right state velocity.
   * @param PR Right state pressure.
   * @param mflux Mass flux solution.
   * @param pflux Momentum flux solution.
   * @param Eflux Energy flux solution.
   * @return Flag signaling whether the left state (-1), the right state (1), or
   * a vacuum state (0) was sampled.
   */
  inline int solve_for_flux(double rhoL, double uL, double PL, double rhoR,
                            double uR, double PR, double &mflux, double &pflux,
                            double &Eflux) {

    // Handle vacuum
    if (rhoL == 0. && rhoR == 0.) {
      mflux = 0.;
      pflux = 0.;
      Eflux = 0.;
      return 0;
    }

    // precompute inverse densities, as they are used multiple times below
    // (and multiplication is much faster than division)
    const double rhoLinv = 1. / rhoL;
    const double rhoRinv = 1. / rhoR;
    const double aL = std::sqrt(_gamma * PL * rhoLinv);
    const double aR = std::sqrt(_gamma * PR * rhoRinv);

    // precompute the velocity difference, which is used multiple times below
    const double uRmuL = uR - uL;

    // Handle vacuum: vacuum does not require iteration and is always exact
    // We haven't implemented this, but should it ever be necessary, it is
    // straightforward to copy this over from the exact Riemann solver.
    if (rhoL == 0. || rhoR == 0.) {
      std::cerr << "Not supported yet!" << std::endl;
      abort();
      return 0;
    }
    if (2. * _odgm1 * (aL + aR) <= uRmuL) {
      std::cerr << "Not supported yet!" << std::endl;
      abort();
      return 0;
    }

    // STEP 1: pressure estimate
    const double Ppvrs =
        0.5 * (PL + PR) - 0.125 * uRmuL * (rhoL + rhoR) * (aL + aR);
    const double Pstar = std::max(0., Ppvrs);

    // STEP 2: wave speed estimates
    // all these speeds are along the interface normal, since uL and uR are
    double qL = 1.;
    if (Pstar > PL) {
      qL = std::sqrt(1. + _hgp1dg * (Pstar / PL - 1.));
    }
    double qR = 1.;
    if (Pstar > PR) {
      qR = std::sqrt(1. + _hgp1dg * (Pstar / PR - 1.));
    }
    const double SL = uL - aL * qL;
    const double SR = uR + aR * qR;
    const double SLmuL = SL - uL;
    const double SRmuR = SR - uR;
    const double Sstar = (PR - PL + rhoL * uL * SLmuL - rhoR * uR * SRmuR) /
                         (rhoL * SLmuL - rhoR * SRmuR);

    if (Sstar >= 0.) {
      // flux FL
      const double uL2 = uL * uL;
      const double rhoLuL = rhoL * uL;
      const double eL = _odgm1 * PL * rhoLinv + 0.5 * uL2;
      const double rhoLeL = rhoL * eL;
      mflux = rhoLuL;
      pflux = rhoL * uL2 + PL;
      Eflux = (rhoLeL + PL) * uL;
      if (SL < 0.) {
        // flux FL*
        const double rhostarL = rhoL * SLmuL / (SL - Sstar);
        const double ustarL = rhostarL * Sstar;
        const double PstarL =
            rhostarL * (eL + (Sstar - uL) * (Sstar + PL / (rhoL * SLmuL)));
        mflux += SL * (rhostarL - rhoL);
        pflux += SL * (ustarL - rhoLuL);
        Eflux += SL * (PstarL - rhoLeL);
      }
      return -1;
    } else {
      // flux FR
      const double uR2 = uR * uR;
      const double rhoRuR = rhoR * uR;
      const double eR = _odgm1 * PR * rhoRinv + 0.5 * uR2;
      const double rhoReR = rhoR * eR;
      mflux = rhoRuR;
      pflux = rhoR * uR2 + PR;
      Eflux = (rhoReR + PR) * uR;
      if (SR > 0.) {
        // flux FR*
        const double rhostarR = rhoR * SRmuR / (SR - Sstar);
        const double ustarR = rhostarR * Sstar;
        const double PstarR =
            rhostarR * (eR + (Sstar - uR) * (Sstar + PR / (rhoR * SRmuR)));
        mflux += SR * (rhostarR - rhoR);
        pflux += SR * (ustarR - rhoRuR);
        Eflux += SR * (PstarR - rhoReR);
      }
      return 1;
    }
  }
};

#endif // HLLCRIEMANNSOLVER_HPP
