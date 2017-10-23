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

#include "Bondi.hpp"
#include "Boundaries.hpp"
#include "Cell.hpp"
#include "EOS.hpp"
#include "IC.hpp"
#include "Potential.hpp"
#include "RiemannSolver.hpp"
#include "SafeParameters.hpp"
#include "Spherical.hpp"
#include "Units.hpp"

// standard libraries
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>

/**
 * @brief Get the neutral fraction for the cell with the given midpoint radius.
 *
 * Old version that contains a strong jump from ionized to neutral.
 *
 * @param rmin Radius of the lower wall of the cell (in internal units of L).
 * @param rmax Radius of the upper wall of the cell (in internal units of L).
 * @param rion Ionization radius (in internal units of L).
 * @param transition_width Width of the transition region (in internal units of
 * L).
 * @return Neutral fraction within the cell.
 */
double get_neutral_fraction2(double rmin, double rmax, double rion,
                             double transition_width) {
  if (rion > rmax) {
    return 0.;
  } else {
    if (rion < rmin) {
      return 1.;
    } else {
      return (rion - rmax) / (rmin - rmax);
    }
  }
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
double get_neutral_fraction3(double rmin, double rmax, double rion,
                             double transition_width) {
  const double rion_min = rion - 0.5 * transition_width;
  const double rion_max = rion + 0.5 * transition_width;
  double slope;
  if (rion_max - rion_min > 0.) {
    slope = 1. / (rion_max - rion_min);
  } else {
    slope = 0.;
  }
  // Note: we assume rmax - rmin << rion_max - rion_min
  if (rmax < rion_min) {
    return 0.;
  } else if (rmin < rion_min && rmax > rion_min) {
    return 0.5 * (rmax - rion_min) * (rmax - rion_min) * slope / (rmax - rmin);
  } else if (rmin > rion_min && rmax < rion_max) {
    return ((rmin - rion_min) * slope * (rmax - rmin) +
            0.5 * (rmax - rmin) * slope * (rmax - rmin)) /
           (rmax - rmin);
  } else if (rmin < rion_max && rmax > rion_max) {
    return ((rmax - rion_max) + (rmin - rion_min) * slope * (rion_max - rmin) +
            0.5 * (rion_max - rmin) * (rion_max - rmin) * slope) /
           (rmax - rmin);
  } else {
    return 1.;
  }
}

double get_neutral_fraction_integral(double A, double S, double rion,
                                     double rmin, double rmax) {
  const double rmaxrel = rmax - rion;
  const double rminrel = rmin - rion;
  const double rdiff = rmax - rmin;
  const double rmaxrel2 = rmaxrel * rmaxrel;
  const double rmaxrel4 = rmaxrel2 * rmaxrel2;
  const double rminrel2 = rminrel * rminrel;
  const double rminrel4 = rminrel2 * rminrel2;
  return 0.25 * A * (rmaxrel4 - rminrel4) + 0.5 * S * (rmaxrel2 - rminrel2) +
         0.5 * rdiff;
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
double get_neutral_fraction(double rmin, double rmax, double rion,
                            double transition_width) {
  const double rion_min = rion - 0.5 * transition_width;
  const double rion_max = rion + 0.5 * transition_width;
  double S, A;
  if (transition_width > 0.) {
    S = 3. / (2. * transition_width);
    A = -16. * S * S * S / 27.;
  } else {
    S = 0.;
    A = 0.;
  }
  // Note: we assume rmax - rmin << rion_max - rion_min
  if (rmax < rion_min) {
    return 0.;
  } else if (rmin < rion_min && rmax > rion_min) {
    return get_neutral_fraction_integral(A, S, rion, rion_min, rmax) /
           (rmax - rmin);
  } else if (rmin > rion_min && rmax < rion_max) {
    return get_neutral_fraction_integral(A, S, rion, rmin, rmax) /
           (rmax - rmin);
  } else if (rmin < rion_max && rmax > rion_max) {
    return (get_neutral_fraction_integral(A, S, rion, rmin, rion_max) +
            (rmax - rion_max)) /
           (rmax - rmin);
  } else {
    return 1.;
  }
}

/**
 * @brief Write a snapshot with the given index.
 *
 * @param istep Index of the snapshot file.
 * @param cells Cells to write.
 */
void write_snapshot(unsigned int istep, const Cell *cells,
                    const unsigned int ncell) {
  std::stringstream filename;
  filename << "snapshot_";
  filename.fill('0');
  filename.width(4);
  filename << (istep / SNAPSTEP);
  filename << ".txt";
  std::ofstream ofile(filename.str().c_str());
  ofile << "# time: " << istep * DT * UNIT_TIME_IN_SI << "\n";
  for (unsigned int i = 1; i < ncell + 1; ++i) {
    ofile << cells[i]._midpoint * UNIT_LENGTH_IN_SI << "\t"
          << cells[i]._rho * UNIT_DENSITY_IN_SI << "\t"
          << cells[i]._u * UNIT_VELOCITY_IN_SI << "\t"
          << cells[i]._P * UNIT_PRESSURE_IN_SI << "\t" << cells[i]._nfac
          << "\n";
  }
  ofile.close();
}

/**
 * @brief Write a binary snapshot that can be used as initial condition file.
 *
 * @param Cells to write.
 */
void write_binary_snapshot(const Cell *cells, const unsigned int ncell) {
  std::ofstream ofile("lastsnap.dat");
  for (unsigned int i = 1; i < ncell + 1; ++i) {
    ofile.write(reinterpret_cast<const char *>(&cells[i]._rho), sizeof(double));
    ofile.write(reinterpret_cast<const char *>(&cells[i]._u), sizeof(double));
    ofile.write(reinterpret_cast<const char *>(&cells[i]._P), sizeof(double));
    ofile.write(reinterpret_cast<const char *>(&cells[i]._a), sizeof(double));
  }
}

/**
 * @brief Main simulation program.
 *
 * @param argc Number of command line arguments (ignored).
 * @param argv Command line arguments (ignored).
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  unsigned int ncell = NCELL;
  std::string ic_file_name(IC_FILE_NAME);
  double transition_width = IONIZATION_TRANSITION_WIDTH;

  if (argc > 1) {
    ncell = atoi(argv[1]);
  }
  if (argc > 2) {
    ic_file_name = argv[2];
  }
  if (argc > 3) {
    transition_width = atof(argv[3]);
  }

  //  for(unsigned int i = 0; i < 1000; ++i){
  //    const double x = 0.2 + 0.001 * 0.2 * i;
  //    std::cout << x << "\t" << get_neutral_fraction(x, 0.3) << std::endl;
  //  }
  //  return 0;
  std::cout << "Slope: " << (1.5 / transition_width) / UNIT_LENGTH_IN_SI
            << std::endl;

  std::cout << "UNIT_LENGTH_IN_SI: " << UNIT_LENGTH_IN_SI << std::endl;
  std::cout << "UNIT_MASS_IN_SI: " << UNIT_MASS_IN_SI << std::endl;
  std::cout << "UNIT_TIME_IN_SI: " << UNIT_TIME_IN_SI << std::endl;
  std::cout << "UNIT_DENSITY_IN_SI: " << UNIT_DENSITY_IN_SI << std::endl;
  std::cout << "UNIT_VELOCITY_IN_SI: " << UNIT_VELOCITY_IN_SI << std::endl;
  std::cout << "UNIT_PRESSURE_IN_SI: " << UNIT_PRESSURE_IN_SI << std::endl;

#if EOS == EOS_ISOTHERMAL || EOS == EOS_BONDI
  std::cout << "Newton G: "
            << G * (UNIT_LENGTH_IN_SI * UNIT_LENGTH_IN_SI * UNIT_LENGTH_IN_SI /
                    UNIT_MASS_IN_SI / UNIT_TIME_IN_SI / UNIT_TIME_IN_SI)
            << " m^3 kg^-1 s^-2" << std::endl;
  std::cout << "ISOTHERMAL_C_SQUARED: " << ISOTHERMAL_C_SQUARED << std::endl;
  std::cout << "Neutral sound speed: "
            << std::sqrt(ISOTHERMAL_C_SQUARED) * UNIT_VELOCITY_IN_SI
            << " m s^-1" << std::endl;
  std::cout << "Neutral temperature: "
            << ISOTHERMAL_C_SQUARED * HYDROGEN_MASS_IN_SI *
                   UNIT_VELOCITY_IN_SI * UNIT_VELOCITY_IN_SI / BOLTZMANN_K_IN_SI
            << " K" << std::endl;
#endif

  std::cout << "Point mass: " << MASS_POINT_MASS * UNIT_MASS_IN_SI << " kg"
            << std::endl;

  std::cout << "Useful units:" << std::endl;
  std::cout << "Point mass: " << MASS_POINT_MASS * UNIT_MASS_IN_MSOL << " Msol"
            << std::endl;
  std::cout << "Time step: " << DT * UNIT_TIME_IN_YR << " yr" << std::endl;
  std::cout << "Total simulation time: " << DT * NSTEP * UNIT_TIME_IN_YR
            << " yr" << std::endl;
  std::cout << "Time in between snapshots: " << DT * SNAPSTEP * UNIT_TIME_IN_YR
            << " yr " << std::endl;
  std::cout << "Minimum radius: " << RMIN * UNIT_LENGTH_IN_AU << " AU"
            << std::endl;
  std::cout << "Maximum radius: " << RMAX * UNIT_LENGTH_IN_AU << " AU"
            << std::endl;

  // create the 1D spherical grid
  // we create 2 ghost cells to the left and to the right of the simulation box
  // to handle boundary conditions
  Cell *cells = new Cell[ncell + 2];
  for (unsigned int i = 0; i < ncell + 2; ++i) {
    cells[i]._midpoint = RMIN + (i - 0.5) * BOXSIZE / ncell;
    cells[i]._V = BOXSIZE / ncell;
    const double xmin = RMIN + (i - 1.) * BOXSIZE / ncell;
    const double xmax = RMIN + i * BOXSIZE / ncell;
    cells[i]._Vreal =
        4. * M_PI * (xmax * xmax * xmax - xmin * xmin * xmin) / 3.;
    if (i < ncell + 1) {
      cells[i]._right_ngb = &cells[i + 1];
    }
  }

  // set up the initial condition
  initialize(cells, ncell);

  // convert primitive variables to conserved variables
  for (unsigned int i = 1; i < ncell + 1; ++i) {
    // apply the equation of state to get the initial pressure (if necessary)
    initial_pressure(cells[i]);

    // use the cell volume to convert primitive into conserved variables
    cells[i]._m = cells[i]._rho * cells[i]._V;
    cells[i]._p = cells[i]._m * cells[i]._u;
    cells[i]._E = cells[i]._P * cells[i]._V / (GAMMA - 1.) +
                  0.5 * cells[i]._u * cells[i]._p;
  }

  // initialize the exact Riemann solver
  RiemannSolver solver(GAMMA);

  // main simulation loop: perform NSTEP steps
  for (unsigned int istep = 0; istep < NSTEP; ++istep) {

    add_spherical_source_term();

    // do first gravity kick
    do_gravity();

    before_primitive_variable_conversion();
    // update the primitive variables based on the values of the conserved
    // variables and the current cell volume
    for (unsigned int i = 1; i < ncell + 1; ++i) {
      cells[i]._rho = cells[i]._m / cells[i]._V;
      cells[i]._u = cells[i]._p / cells[i]._m;
      update_pressure(cells[i]);
    }
    after_primitive_variable_conversion();

    double maxcs = 0.;
    for (unsigned int i = 1; i < ncell + 1; ++i) {
      const double cs = std::sqrt(GAMMA * cells[i]._P / cells[i]._rho) +
                        std::abs(cells[i]._u);
      maxcs = std::max(maxcs, cs);
    }
    const double maxdt = 0.5 * CELLSIZE / maxcs;
    if (maxdt < DT) {
      std::cout << "Time step too large: " << maxdt << " < " << DT << "!"
                << std::endl;
      return 1;
    }

    // check if we need to output a snapshot
    if (istep % SNAPSTEP == 0) {
      const double pct = istep * 100. / NSTEP;
      std::cout << "step " << istep << " of " << NSTEP << " (" << pct << " %)"
                << std::endl;
      std::cout << "Maximal system time step: " << maxdt << std::endl;
      write_snapshot(istep, cells, ncell);
    }

    // apply boundary conditions
    boundary_conditions_primitive_variables();

    // compute slope limited gradients for the primitive variables in each cell
    for (unsigned int i = 1; i < ncell + 1; ++i) {
      const double dx = cells[i + 1]._midpoint - cells[i - 1]._midpoint;
      const double dx_inv = 1. / dx;

      const double gradrho = (cells[i + 1]._rho - cells[i - 1]._rho) * dx_inv;
      const double rhomax = std::max(cells[i - 1]._rho, cells[i + 1]._rho);
      const double rhomin = std::min(cells[i - 1]._rho, cells[i + 1]._rho);
      const double rho_ext_plu = HALF_CELLSIZE * gradrho;
      const double rho_ext_min = -HALF_CELLSIZE * gradrho;
      const double rhoextmax = std::max(rho_ext_min, rho_ext_plu);
      const double rhoextmin = std::min(rho_ext_min, rho_ext_plu);
      const double alpha_rho =
          std::min(1., 0.5 * std::min((rhomax - cells[i]._rho) / rhoextmax,
                                      (rhomin - cells[i]._rho) / rhoextmin));
      cells[i]._grad_rho = alpha_rho * gradrho;

      const double gradu = (cells[i + 1]._u - cells[i - 1]._u) * dx_inv;
      const double umax = std::max(cells[i - 1]._u, cells[i + 1]._u);
      const double umin = std::min(cells[i - 1]._u, cells[i + 1]._u);
      const double u_ext_plu = HALF_CELLSIZE * gradu;
      const double u_ext_min = -HALF_CELLSIZE * gradu;
      const double uextmax = std::max(u_ext_min, u_ext_plu);
      const double uextmin = std::min(u_ext_min, u_ext_plu);
      const double alpha_u =
          std::min(1., 0.5 * std::min((umax - cells[i]._u) / uextmax,
                                      (umin - cells[i]._u) / uextmin));
      cells[i]._grad_u = alpha_u * gradu;

      const double gradP = (cells[i + 1]._P - cells[i - 1]._P) * dx_inv;
      const double Pmax = std::max(cells[i - 1]._P, cells[i + 1]._P);
      const double Pmin = std::min(cells[i - 1]._P, cells[i + 1]._P);
      const double P_ext_plu = HALF_CELLSIZE * gradP;
      const double P_ext_min = -HALF_CELLSIZE * gradP;
      const double Pextmax = std::max(P_ext_min, P_ext_plu);
      const double Pextmin = std::min(P_ext_min, P_ext_plu);
      const double alpha_P =
          std::min(1., 0.5 * std::min((Pmax - cells[i]._P) / Pextmax,
                                      (Pmin - cells[i]._P) / Pextmin));
      cells[i]._grad_P = alpha_P * gradP;
    }

    // apply boundary conditions
    boundary_conditions_gradients();

#ifdef NO_GRADIENTS
    // reset all gradients to zero to disable the second order scheme
    for (unsigned int i = 0; i < ncell + 2; ++i) {
      cells[i]._grad_rho = 0.;
      cells[i]._grad_u = 0.;
      cells[i]._grad_P = 0.;
    }
#endif

    // evolve all primitive variables forward in time for half a time step
    // using the Euler equations and the spatial gradients within the cells
    for (unsigned int i = 0; i < ncell + 2; ++i) {
      const double rho = cells[i]._rho;
      const double u = cells[i]._u;
      const double P = cells[i]._P;
      cells[i]._rho -=
          0.5 * DT * (rho * cells[i]._grad_u + u * cells[i]._grad_rho);
      cells[i]._u -= 0.5 * DT * (u * cells[i]._grad_u + cells[i]._grad_P / rho);
      cells[i]._P -=
          0.5 * DT * (GAMMA * P * cells[i]._grad_u + u * cells[i]._grad_P);

      if (cells[i]._rho < 0.) {
        cells[i]._rho = rho;
      }
      if (cells[i]._P < 0.) {
        cells[i]._P = P;
      }

      add_gravitational_prediction(cells[i]);
    }

    // do the flux exchange
    // we always do the flux for the right face of the cell, which means we
    // start from the left ghost and do not include the right ghost
    for (unsigned int i = 0; i < ncell + 1; ++i) {
      // get the variables in the left and right state
      double rhoL, uL, PL, rhoR, uR, PR;
      rhoL = cells[i]._rho;
      uL = cells[i]._u;
      PL = cells[i]._P;
      Cell *rcell = cells[i]._right_ngb;
      rhoR = rcell->_rho;
      uR = rcell->_u;
      PR = rcell->_P;

      // do the second order spatial reconstruction
      double dmin, dplu;
      dmin = 0.5 * (rcell->_midpoint - cells[i]._midpoint);
      dplu = -dmin;
      double rhoL_dash, rhoR_dash, uL_dash, uR_dash, PL_dash, PR_dash;
      rhoL_dash = rhoL + dmin * cells[i]._grad_rho;
      uL_dash = uL + dmin * cells[i]._grad_u;
      PL_dash = PL + dmin * cells[i]._grad_P;
      rhoR_dash = rhoR + dplu * rcell->_grad_rho;
      uR_dash = uR + dplu * rcell->_grad_u;
      PR_dash = PR + dplu * rcell->_grad_P;

      if (rhoL_dash < 0.) {
        rhoL_dash = rhoL;
      }
      if (rhoR_dash < 0.) {
        rhoR_dash = rhoR;
      }
      if (PL_dash < 0.) {
        PL_dash = PL;
      }
      if (PR_dash < 0.) {
        PR_dash = PR;
      }

      // solve the Riemann problem across the cell face
      double rhosol, usol, Psol;
      solver.solve(rhoL_dash, uL_dash, PL_dash, rhoR_dash, uR_dash, PR_dash,
                   rhosol, usol, Psol);

      // compute and exchange fluxes
      double mflux, pflux, Eflux;
      mflux = rhosol * usol;
      pflux = rhosol * usol * usol + Psol;
      Eflux = (Psol / (GAMMA - 1.) + 0.5 * rhosol * usol * usol) * usol +
              Psol * usol;
      cells[i]._m -= DT * mflux;
      cells[i]._p -= DT * pflux;
      cells[i]._E -= DT * Eflux;
      rcell->_m += DT * mflux;
      rcell->_p += DT * pflux;
      rcell->_E += DT * Eflux;
    }

    add_spherical_source_term();

    // do second gravity kick
    do_gravity();
  }

  // write the final snapshot
  write_snapshot(NSTEP, cells, ncell);
  write_binary_snapshot(cells, ncell);

  delete[] cells;

  // all went well: return with exit code 0
  return 0;
}
