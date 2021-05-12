/*******************************************************************************
 * This file is part of HydroCodeSpherical1D
 * Copyright (C) 2021 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file Cooling.hpp
 *
 * @brief Cooling function.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef COOLING_HPP
#define COOLING_HPP

#define SPECIFIC_ENERGY_CONVERSION_FACTOR                                      \
  ((UNIT_LENGTH_IN_SI * UNIT_LENGTH_IN_SI * UNIT_LENGTH_IN_SI) /               \
   UNIT_ENERGY_IN_SI)

#if COOLING == COOLING_CURVE

#if DIMENSIONALITY == DIMENSIONALITY_1D
#error "Cooling from a curve only works in 3D!"
#endif

// Cooling table: hard-coded, since it is small and then we don't need to read
// files and make sure they can be found.

/*! @brief Cooling table from Kirchschlager (2019): temperatures (in K). */
const double cooling_T[73] = {10000.0,
                              10290.401663604674,
                              11213.193261435596,
                              12102.697512665345,
                              13898.622482052577,
                              14694.1925858846,
                              16062.904513624355,
                              17392.35209983816,
                              18952.01098477871,
                              21049.437649811163,
                              24057.91704843993,
                              28026.170919611006,
                              32962.00480472719,
                              38034.28465959962,
                              43292.90526959403,
                              47368.51374500882,
                              52111.17284883071,
                              60130.176966474435,
                              71853.70760119933,
                              88356.44005686993,
                              111804.57203085616,
                              135313.69158224564,
                              171223.3920643563,
                              206699.85868755623,
                              268978.51195546647,
                              309384.00077478465,
                              336592.2769152817,
                              354729.0420671322,
                              377787.7444391393,
                              397764.5304142671,
                              422006.87766226433,
                              450297.4315273315,
                              472302.24999443145,
                              509768.37132494425,
                              550206.5518560404,
                              588212.8669143564,
                              634873.7808529241,
                              698438.9991485296,
                              772043.2411130213,
                              909870.3245098612,
                              1161843.4799728123,
                              1362122.3143026093,
                              1542018.2863407205,
                              1830979.512961513,
                              2132992.0074644615,
                              2631234.9245853107,
                              3380200.680058953,
                              3966032.873523129,
                              4620213.585505019,
                              5637626.827085784,
                              6897864.017584717,
                              8649464.333106754,
                              10988720.955436977,
                              13744422.600622471,
                              16549642.645536577,
                              19704554.153339442,
                              27781426.05193267,
                              37433284.87932664,
                              46508937.47758922,
                              56288360.19267283,
                              74944862.17586084,
                              92451138.15957102,
                              114046683.10609132,
                              140686704.20316154,
                              173549534.28269237,
                              214088751.45901987,
                              264097473.32783756,
                              325787669.56609106,
                              401888001.0622067,
                              495764513.16556245,
                              611569521.5201099,
                              754425275.9524019,
                              887291012.4639999};

/*! @brief Cooling table from Kirchschlager (2019): cooling function (originally
 *  expressed in erg cm^3 s^-1, converted to J m^3 s^-1). */
const double cooling_Lambda[73] = {0.0,
                                   9.89256025093041e-35,
                                   1.2464089281312487e-34,
                                   1.5945833729313225e-34,
                                   2.6934252937966626e-34,
                                   3.5561163107809447e-34,
                                   4.699608244728035e-34,
                                   6.879819455701967e-34,
                                   9.162286199293358e-34,
                                   1.1597013176960812e-33,
                                   1.493604350984838e-33,
                                   1.9163122741823547e-33,
                                   2.51004589290006e-33,
                                   3.233771336327001e-33,
                                   4.197351912521545e-33,
                                   5.172856347506641e-33,
                                   6.645714455511262e-33,
                                   9.668161329962818e-33,
                                   1.4815098119339133e-32,
                                   2.1914023573069507e-32,
                                   3.293449443129264e-32,
                                   4.5276900190861674e-32,
                                   6.098961992770102e-32,
                                   7.209536533426396e-32,
                                   6.957950911672741e-32,
                                   4.4590449779941983e-32,
                                   2.7269876428756634e-32,
                                   2.010949616091764e-32,
                                   1.3889044303750353e-32,
                                   9.901723964527769e-33,
                                   7.206199669019085e-33,
                                   5.133493324284437e-33,
                                   4.113466122740202e-33,
                                   2.944902738259438e-33,
                                   2.1123376546894048e-33,
                                   1.6198810306582637e-33,
                                   1.200716269830946e-33,
                                   8.873004556347361e-34,
                                   6.872316061811136e-34,
                                   4.982143094279962e-34,
                                   4.910625351114225e-34,
                                   6.1635360506311285e-34,
                                   8.207918904296672e-34,
                                   1.1205152631170961e-33,
                                   1.3943725423062007e-33,
                                   1.4161086984020296e-33,
                                   1.0930402967726974e-33,
                                   8.731807042763502e-34,
                                   6.678252157185632e-34,
                                   4.93345245600413e-34,
                                   3.775394909046169e-34,
                                   2.9076246242681243e-34,
                                   2.348768666513851e-34,
                                   2.0312709830969887e-34,
                                   1.8236394446008286e-34,
                                   1.5945833729313225e-34,
                                   1.4660731852650456e-34,
                                   1.3660592880038619e-34,
                                   1.3192131964215211e-34,
                                   1.3040330499558834e-34,
                                   1.3192131964215211e-34,
                                   1.3554288841281864e-34,
                                   1.3427815380314463e-34,
                                   1.4249208283704266e-34,
                                   1.5231502029381332e-34,
                                   1.6469137494602603e-34,
                                   1.784447305845701e-34,
                                   1.9455878715999643e-34,
                                   2.1323566629225293e-34,
                                   2.346812239959394e-34,
                                   2.590919780479011e-34,
                                   2.8278310840725746e-34,
                                   3.0933788164782573e-34};

/**
 * @brief Get the cooling for the given density and specific thermal energy.
 *
 * @param rho Density.
 * @param etherm Specific thermal energy.
 * @return Change of specific thermal energy over time.
 */
static inline double get_cooling(const double rho, const double etherm) {

  // compute the number density and temperature
  const double n = rho / HYDROGEN_MASS_IN_SI;
  const double T = (GAMMA - 1.) * etherm / BOLTZMANN_K_IN_SI / n;

  // now figure out the index of the highest tabulated temperature below the
  // cell temperature
  uint_fast32_t jl;
  if (T <= cooling_T[0]) {
    // temperature is too low: no cooling
    return 0.;
  } else if (T >= cooling_T[72]) {
    // temperature is too high: linearly extrapolate the tail of the table
    jl = 71;
  } else {
    // temperature within range: use bisection to find the right bin
    jl = 0;
    uint_fast32_t ju = 73;
    while (ju - jl > 1) {
      const uint_fast32_t jm = (ju + jl) >> 1;
      if (T > cooling_T[jm]) {
        jl = jm;
      } else {
        ju = jm;
      }
    }
    // out of range (should have been caught before): revert to extrapolation
    if (jl == 72) {
      jl = 71;
    }
  }

  // now linearly inter-/extrapolate to find the right cooling
  const double Tm = cooling_T[jl];
  const double Tp = cooling_T[jl + 1];
  const double Lambdam = cooling_Lambda[jl];
  const double Lambdap = cooling_Lambda[jl + 1];
  const double fT = (T - Tm) / (Tp - Tm);
  const double Lambda = Lambdam * (1. - fT) + Lambdap * fT;

  // normalise with the cell density squared (and apply the right sign)
  return -Lambda * n * n;
}

/**
 * @brief Explicitly compute the new specific thermal energy due to cooling over
 * the given time step.
 *
 * This function a fourth order Runge-Kutta scheme, Eq (15.37) in Toro (2009).
 *
 * @param dt Time step.
 * @param rho Density.
 * @param etherm Specific thermal energy.
 * @return New specific thermal energy.
 */
static inline double cooling_fourth_order_runge_kutta(const double dt,
                                                      const double rho,
                                                      const double etherm) {
  const double K1 = dt * get_cooling(rho, etherm);
  const double K2 = dt * get_cooling(rho, etherm + 0.5 * K1);
  const double K3 = dt * get_cooling(rho, etherm + 0.5 * K2);
  const double K4 = dt * get_cooling(rho, etherm + K3);
  double etherm_new = etherm + (K1 + 2. * K2 + 2. * K3 + K4) / 6.;
  return etherm_new;
}

/**
 * @brief Compute the new energy after cooling over the given timestep.
 *
 * @param dt Time step.
 * @param U Specific conserved variables.
 * @return New value for the specific conserved energy.
 */
static inline double energy_after_cooling(const double dt, const double U[3]) {
  // precompute the specific kinetic energy
  const double ekin = 0.5 * U[1] * U[1] / U[0];
  // get the specific thermal energy
  double etherm = U[2] - ekin;
  // convert to SI units, since the table uses these
  etherm /= SPECIFIC_ENERGY_CONVERSION_FACTOR;
  // compute the new specific thermal energy
  etherm = cooling_fourth_order_runge_kutta(dt * UNIT_TIME_IN_SI,
                                            U[0] * UNIT_DENSITY_IN_SI, etherm);
  // convert back to code units
  etherm *= SPECIFIC_ENERGY_CONVERSION_FACTOR;
  // apply the minimum temperature limit
  etherm = std::max(etherm, U[0] * MINIMUM_SPECIFIC_THERMAL_ENERGY);
  // add back the kinetic energy and return
  return etherm + ekin;
}

#endif // COOLING == COOLING_CURVE

/**
 * @brief Update the energy with the cooling contribution.
 */
#if COOLING == COOLING_NONE
#define add_cooling()
#elif COOLING == COOLING_CURVE
#define add_cooling()                                                          \
  _Pragma("omp parallel for") for (uint_fast32_t i = 1; i < ncell + 1; ++i) {  \
    if (cells[i]._m > 0.) {                                                    \
      const double Vinv = 1. / cells[i]._V;                                    \
      const double dt = cells[i]._dt;                                          \
      const double Ui[3] = {cells[i]._m * Vinv, cells[i]._p * Vinv,            \
                            cells[i]._E * Vinv};                               \
      cells[i]._E = energy_after_cooling(0.5 * dt, Ui) * cells[i]._V;          \
    }                                                                          \
  }
#endif

#endif // COOLING_HPP
