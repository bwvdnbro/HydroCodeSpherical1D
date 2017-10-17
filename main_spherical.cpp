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

// exact Riemann solver
#include "RiemannSolver.hpp"

// standard libraries
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>

// USEFUL SIMULATION CONSTANTS
// These are used as aliases to check in the rest of the program.
// Since invalid defines default to 0, we have to start from 1 (otherwise our
// sensibility checks don't work)

// Possible types of boundary conditions.
/*! @brief Open boundaries: inflow or outflow depending on the local flow
 *  velocity */
#define BOUNDARIES_OPEN 1
/*! @brief Reflective boundaries: outgoing flows are reflected inwards. */
#define BOUNDARIES_REFLECTIVE 2
/*! @brief Bondi boundaries: impose the Bondi solution on the boundaries. */
#define BOUNDARIES_BONDI 3

// Possible types of equation of state.
/*! @brief Ideal gas (adiabatic) equation of state. */
#define EOS_IDEAL 1
/*! @brief Isothermal equation of state. */
#define EOS_ISOTHERMAL 2
/*! @brief Bondi equation of state: isothermal gas with ionization pressure. */
#define EOS_BONDI 3

// Possible types of external potentials.
/*! @brief No external potential (no gravity). */
#define POTENTIAL_NONE 1
/*! @brief Point mass external potential. */
#define POTENTIAL_POINT_MASS 2

// Possible types of initial conditions
/*! @brief 1D spherical Sod shock setup. */
#define IC_SOD 1
/*! @brief 1D spherical Bondi accretion setup. */
#define IC_BONDI 2

// PHYSICAL CONSTANTS

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

// UNIT CONVERSIONS

/*! @brief Time unit (in s). */
#define UNIT_TIME_IN_SI                                                        \
  std::sqrt(G *UNIT_LENGTH_IN_SI *UNIT_LENGTH_IN_SI *UNIT_LENGTH_IN_SI /       \
            (UNIT_MASS_IN_SI * NEWTON_G_IN_SI))

/*! @brief Density unit (in kg m^-3). */
#define UNIT_DENSITY_IN_SI                                                     \
  (UNIT_MASS_IN_SI /                                                           \
   (UNIT_LENGTH_IN_SI * UNIT_LENGTH_IN_SI * UNIT_LENGTH_IN_SI))

/*! @brief Velocity unit (in m s^-1). */
#define UNIT_VELOCITY_IN_SI (UNIT_LENGTH_IN_SI / UNIT_TIME_IN_SI)

/*! @brief Pressure unit (in kg m^-2 s^-2). */
#define UNIT_PRESSURE_IN_SI                                                    \
  (UNIT_MASS_IN_SI / (UNIT_LENGTH_IN_SI * UNIT_LENGTH_IN_SI *                  \
                      UNIT_TIME_IN_SI * UNIT_TIME_IN_SI))

/*! @brief Mass unit (in Msol). */
#define UNIT_MASS_IN_MSOL (UNIT_MASS_IN_SI / SOLAR_MASS_IN_SI)

/*! @brief Time unit (in yr). */
#define UNIT_TIME_IN_YR (UNIT_TIME_IN_SI / YEAR_IN_SI)

/*! @brief Lenght unit (in AU). */
#define UNIT_LENGTH_IN_AU (UNIT_LENGTH_IN_SI / AU_IN_SI)

// READ SIMULATION PARAMETERS
// for clarity, they are stored in a separate file
#include "Parameters.hpp"

// end of customizable parameters, the parameters below depend on the ones
// above

/*! @brief Size of the simulation "box" (in internal units of L). */
#define BOXSIZE (RMAX - RMIN)
/*! @brief Half the size of a single "cell" of the simulation (in internal units
 *  of L). */
#define CELLSIZE (BOXSIZE / NCELL)
#define HALF_CELLSIZE (0.5 * CELLSIZE)
/*! @brief Isothermal sound speed squared (if EOS_ISOTHERMAL is selected, in
 *  internal units of L T^-1). */
#define ISOTHERMAL_C_SQUARED                                                   \
  (ISOTHERMAL_TEMPERATURE * BOLTZMANN_K_IN_SI / HYDROGEN_MASS_IN_SI /          \
   UNIT_VELOCITY_IN_SI / UNIT_VELOCITY_IN_SI)

// sanity checks on selected types

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

/**
 * @brief Single cell of the grid.
 */
class Cell {
public:
  // primitive (volume dependent) variables

  /*! @brief Density (in internal units of M L^-3). */
  double _rho;
  /*! @brief Fluid velocity (in internal units of L T^-1). */
  double _u;
  /*! @brief Pressure (in internal units of M L^-1 T^-2). */
  double _P;

  /*! @brief Density gradient (in internal units of M L^-4). */
  double _grad_rho;
  /*! @brief Velocity gradient (in internal units of T^-1). */
  double _grad_u;
  /*! @brief Pressure gradient (in internal units of M L^-2 T^-2). */
  double _grad_P;

  // conserved variables

  /*! @brief Mass (in internal units of M). */
  double _m;
  /*! @brief Momentum (in internal units of L M T^-1). */
  double _p;
  /*! @brief Total energy (in internal units of M L^2 T^-2). */
  double _E;

  // geometrical quantities

  /*! @brief 1D cell volume: radial length of the spherical shell (in internal
   *  units of L). */
  double _V;
  /*! @brief Actual 3D cell volume: volume of the shell (in internal units of
   *  L^3). */
  double _Vreal;
  /*! @brief Radial coordinate of the center of the shell (in internal units of
   *  L). */
  double _midpoint;

  /*! @brief Neighbouring cell to the right. */
  Cell *_right_ngb;

  // gravitational quantities

  /*! @brief Gravitational acceleration (in internal units of L T^-2). */
  double _a;

  // ionization quantities
  double _nfac;
};

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

/**
 * @brief Initialize a cell of the grid.
 *
 * @param cell Cell to initialize.
 */
void init(Cell &cell) {
#if IC == IC_SOD
  if (cell._midpoint < 0.25) {
    cell._rho = 1.;
    cell._P = 1.;
  } else {
    cell._rho = 0.125;
    cell._P = 0.1;
  }
  cell._u = 0.;
  cell._a = 0.;
  cell._nfac = 0.;
#elif IC == IC_BONDI
  const double r_inv = 0.9 / cell._midpoint;
  cell._rho = bondi_density(r_inv);
  cell._u = bondi_velocity(r_inv);
  cell._P = bondi_pressure(r_inv);
  const double r2 = cell._midpoint * cell._midpoint;
  cell._a = -G * MASS_POINT_MASS / r2;
  cell._nfac = 0.;
#endif
}

double get_neutral_fraction2(double r, double rion) {
  const double rmax = r + HALF_CELLSIZE;
  if (rion > rmax) {
    return 0.;
  } else {
    const double rmin = r - HALF_CELLSIZE;
    if (rion < rmin) {
      return 1.;
    } else {
      return (rion - rmax) / (rmin - rmax);
    }
  }
}

double get_neutral_fraction(double r, double rion) {
  const double rmin = r - HALF_CELLSIZE;
  const double rmax = r + HALF_CELLSIZE;
  const double rion_min = rion - 0.5 * IONIZATION_TRANSITION_WIDTH;
  const double rion_max = rion + 0.5 * IONIZATION_TRANSITION_WIDTH;
  const double slope = 1. / (rion_max - rion_min);
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

/**
 * @brief Write a snapshot with the given index.
 *
 * @param istep Index of the snapshot file.
 * @param cells Cells to write.
 */
void write_snapshot(unsigned int istep, const Cell cells[NCELL + 2]) {
  std::stringstream filename;
  filename << "snapshot_";
  filename.fill('0');
  filename.width(4);
  filename << (istep / SNAPSTEP);
  filename << ".txt";
  std::ofstream ofile(filename.str().c_str());
  ofile << "# time: " << istep * DT * UNIT_TIME_IN_SI << "\n";
  for (unsigned int i = 1; i < NCELL + 1; ++i) {
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
void write_binary_snapshot(const Cell cells[NCELL + 2]) {
  std::ofstream ofile("lastsnap.dat");
  for (unsigned int i = 1; i < NCELL + 1; ++i) {
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

  //  for(unsigned int i = 0; i < 1000; ++i){
  //    const double x = 0.2 + 0.001 * 0.2 * i;
  //    std::cout << x << "\t" << get_neutral_fraction(x, 0.3) << std::endl;
  //  }
  //  return 0;

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
  Cell cells[NCELL + 2];
  for (unsigned int i = 0; i < NCELL + 2; ++i) {
    cells[i]._midpoint = RMIN + (i - 0.5) * BOXSIZE / NCELL;
    cells[i]._V = BOXSIZE / NCELL;
    const double xmin = RMIN + (i - 1.) * BOXSIZE / NCELL;
    const double xmax = RMIN + i * BOXSIZE / NCELL;
    cells[i]._Vreal =
        4. * M_PI * (xmax * xmax * xmax - xmin * xmin * xmin) / 3.;
    if (i < NCELL + 1) {
      cells[i]._right_ngb = &cells[i + 1];
    }
  }

#if IC == IC_FILE
  // open the initial condition file
  std::ifstream icfile("ic.dat");
#endif
  // initialize the cells (we don't initialize the ghost cells)
  for (unsigned int i = 1; i < NCELL + 1; ++i) {
#if IC == IC_FILE
    // read the primitive variables from the initial condition file
    icfile.read(reinterpret_cast<char *>(&cells[i]._rho), sizeof(double));
    icfile.read(reinterpret_cast<char *>(&cells[i]._u), sizeof(double));
    icfile.read(reinterpret_cast<char *>(&cells[i]._P), sizeof(double));
    icfile.read(reinterpret_cast<char *>(&cells[i]._a), sizeof(double));
#else
    // get the primitive variables from the initial condition
    init(cells[i]);

#if EOS == EOS_ISOTHERMAL
    // if an isothermal equation of state is used, the pressure is a function of
    // the density, and the initial condition is overwritten
    cells[i]._P = ISOTHERMAL_C_SQUARED * cells[i]._rho;
#endif

#endif

    // use the cell volume to convert primitive into conserved variables
    cells[i]._m = cells[i]._rho * cells[i]._V;
    cells[i]._p = cells[i]._m * cells[i]._u;
    cells[i]._E = cells[i]._P * cells[i]._V / (GAMMA - 1.) +
                  0.5 * cells[i]._u * cells[i]._p;
  }

#if IC == IC_FILE
  // close the initial condition file
  icfile.close();
#endif

  // initialize the exact Riemann solver
  RiemannSolver solver(GAMMA);

  // main simulation loop: perform NSTEP steps
  for (unsigned int istep = 0; istep < NSTEP; ++istep) {

    // add spherical source term (see Toro, 2009, chapter 17)
    // we use a second order Runge-Kutta step and apply an operator splitting
    // method to couple the source term to the hydro step
    for (unsigned int i = 1; i < NCELL + 1; ++i) {
      const double r = cells[i]._midpoint;
      const double Ui[3] = {cells[i]._m / cells[i]._V,
                            cells[i]._p / cells[i]._V,
                            cells[i]._E / cells[i]._V};
      double K1[3];
      K1[0] = -DT * Ui[1] / r;
      K1[1] = -DT * Ui[1] * Ui[1] / Ui[0] / r;
      const double p1 = (GAMMA - 1.) * (Ui[2] - 0.5 * Ui[1] * Ui[1] / Ui[0]);
      K1[2] = -DT * Ui[1] * (Ui[2] + p1) / Ui[0] / r;
      double U[3] = {Ui[0] + K1[0], Ui[1] + K1[1], Ui[2] + K1[2]};
      double K2[3];
      K2[0] = -DT * U[1] / r;
      K2[1] = -DT * U[1] * U[1] / U[0] / r;
      const double p2 = (GAMMA - 1.) * (U[2] - 0.5 * U[1] * U[1] / U[0]);
      K2[2] = -DT * U[1] * (U[2] + p2) / U[0] / r;
      U[0] = Ui[0] + 0.5 * (K1[0] + K2[0]);
      U[1] = Ui[1] + 0.5 * (K1[1] + K2[1]);
      U[2] = Ui[2] + 0.5 * (K1[2] + K2[2]);

      cells[i]._m = U[0] * cells[i]._V;
      cells[i]._p = U[1] * cells[i]._V;
      cells[i]._E = U[2] * cells[i]._V;
    }

#if POTENTIAL == POTENTIAL_POINT_MASS
    // add gravitational acceleration
    for (unsigned int i = 1; i < NCELL + 1; ++i) {
      const double r = cells[i]._midpoint;
      const double a = -G * MASS_POINT_MASS / (r * r);
      cells[i]._a = a;
      const double m = cells[i]._V * cells[i]._rho;
      cells[i]._p += 0.5 * DT * a * m;
      // we do not update the total energy, as we only run gravity simulations
      // with an isothermal eos, in which case the total energy is ignored by
      // the hydro scheme
    }
#endif

#if EOS == EOS_BONDI
    // total ionizing budget of the central source
    // every shell will absorb a fraction of this budget until the ionization
    // radius is reached
    double Cion =
        5. * std::pow(0.9, 2.4) / 3. *
        (std::pow(INITIAL_IONIZATION_RADIUS, 0.6) - std::pow(0.1, 0.6));
    double rion = 0.;
#endif
    // update the primitive variables based on the values of the conserved
    // variables and the current cell volume
    for (unsigned int i = 1; i < NCELL + 1; ++i) {
      cells[i]._rho = cells[i]._m / cells[i]._V;
      cells[i]._u = cells[i]._p / cells[i]._m;
#if EOS == EOS_IDEAL
      cells[i]._P =
          (GAMMA - 1.) * (cells[i]._E / cells[i]._V -
                          0.5 * cells[i]._rho * cells[i]._u * cells[i]._u);
#elif EOS == EOS_ISOTHERMAL
      cells[i]._P = ISOTHERMAL_C_SQUARED * cells[i]._rho;
#elif EOS == EOS_BONDI
      // if ionization is active: check if the cell is ionized
      if (istep > NSTEP_RELAX && Cion >= 0.) {
        // as long as there is ionizing radiation left we assume cells are
        // ionized
        const double rmin = cells[i]._midpoint - HALF_CELLSIZE;
        const double rmax = cells[i]._midpoint + HALF_CELLSIZE;
        const double Vshell = (rmax * rmax * rmax - rmin * rmin * rmin) / 3.;
        const double Cshell = Vshell * cells[i]._rho * cells[i]._rho;
        const double ifac = std::min(1., Cion / Cshell);
        if (ifac < 1.) {
          if (ifac > 0.) {
            const double nfac = 1. - ifac;
            rion = ifac * rmax + nfac * rmin;
          }
        } else {
          rion = rmax;
        }
        // subtract this shell's ionization budget from the total
        Cion -= ifac * Cshell;
      }
#endif
    }
#if EOS == EOS_BONDI
    for (unsigned int i = 1; i < NCELL + 1; ++i) {
      cells[i]._nfac = get_neutral_fraction(cells[i]._midpoint, rion);
    }
    for (unsigned int i = 1; i < NCELL + 1; ++i) {
      const double nfac = cells[i]._nfac;
      const double ifac = 1. - nfac;
      cells[i]._P = ISOTHERMAL_C_SQUARED * cells[i]._rho * (100. * ifac + nfac);
    }
#endif

    // check if we need to output a snapshot
    if (istep % SNAPSTEP == 0) {
      std::cout << "step " << istep << " of " << NSTEP << std::endl;
      write_snapshot(istep, cells);
    }

// apply boundary conditions
#if BOUNDARIES == BOUNDARIES_BONDI
    // impose the Bondi solution at the boundaries
    const double r_inv_low = 0.9 / cells[0]._midpoint;
    cells[0]._rho = bondi_density(r_inv_low);
    cells[0]._u = bondi_velocity(r_inv_low);
    if (istep < NSTEP_RELAX) {
      cells[0]._P = bondi_pressure(r_inv_low);
    } else {
      cells[0]._P = 100. * bondi_pressure(r_inv_low);
    }

    const double r_inv_high = 0.9 / cells[NCELL + 1]._midpoint;
    cells[NCELL + 1]._rho = bondi_density(r_inv_high);
    cells[NCELL + 1]._u = bondi_velocity(r_inv_high);
    cells[NCELL + 1]._P = bondi_pressure(r_inv_high);
#elif BOUNDARIES == BOUNDARIES_OPEN
    cells[0]._rho = cells[1]._rho;
    cells[0]._u = cells[1]._u;
    cells[0]._P = ISOTHERMAL_C_SQUARED * cells[1]._rho;
    cells[NCELL + 1]._rho = cells[NCELL]._rho;
    cells[NCELL + 1]._u = cells[NCELL]._u;
    cells[NCELL + 1]._P = ISOTHERMAL_C_SQUARED * cells[NCELL + 1]._rho;
#elif BOUNDARIES == BOUNDARIES_REFLECTIVE
    cells[0]._rho = cells[1]._rho;
    cells[0]._u = -cells[1]._u;
    cells[0]._P = ISOTHERMAL_C_SQUARED * cells[1]._rho;
    cells[NCELL + 1]._rho = cells[NCELL]._rho;
    cells[NCELL + 1]._u = -cells[NCELL]._u;
    cells[NCELL + 1]._P = cells[NCELL]._P;
#endif

    // compute slope limited gradients for the primitive variables in each cell
    for (unsigned int i = 1; i < NCELL + 1; ++i) {
      double rhomin, rhoplu, umin, uplu, Pmin, Pplu, dmin, dplu;
      rhomin = cells[i]._rho - cells[i - 1]._rho;
      rhoplu = cells[i]._rho - cells[i + 1]._rho;
      umin = cells[i]._u - cells[i - 1]._u;
      uplu = cells[i]._u - cells[i + 1]._u;
      Pmin = cells[i]._P - cells[i - 1]._P;
      Pplu = cells[i]._P - cells[i + 1]._P;
      dmin = cells[i]._midpoint - cells[i - 1]._midpoint;
      dplu = cells[i]._midpoint - cells[i + 1]._midpoint;
      if (std::abs(rhomin * dplu) < std::abs(rhoplu * dmin)) {
        cells[i]._grad_rho = rhomin / dmin;
      } else {
        cells[i]._grad_rho = rhoplu / dplu;
      }
      if (std::abs(umin * dplu) < std::abs(uplu * dmin)) {
        cells[i]._grad_u = umin / dmin;
      } else {
        cells[i]._grad_u = uplu / dplu;
      }
      if (std::abs(Pmin * dplu) < std::abs(Pplu * dmin)) {
        cells[i]._grad_P = Pmin / dmin;
      } else {
        cells[i]._grad_P = Pplu / dplu;
      }
    }

// apply boundary conditions
#if BOUNDARIES == BOUNDARIES_OPEN
    cells[0]._grad_rho = cells[1]._grad_rho;
    cells[0]._grad_u = cells[1]._grad_u;
    cells[0]._grad_P = cells[1]._grad_P;
    cells[NCELL + 1]._grad_rho = cells[NCELL + 1]._grad_rho;
    cells[NCELL + 1]._grad_u = cells[NCELL + 1]._grad_u;
    cells[NCELL + 1]._grad_P = cells[NCELL + 1]._grad_P;
#elif BOUNDARIES == BOUNDARIES_BONDI
    // compute the exact value of the gradient for the bondi solution
    const double rmin = cells[0]._midpoint - CELLSIZE;
    const double rmax = cells[NCELL + 1]._midpoint + CELLSIZE;
    const double rmin_inv = 1. / rmin;
    const double rmax_inv = 1. / rmax;

    // lower boundary
    {
      double rhomin, rhoplu, umin, uplu, Pmin, Pplu, dmin, dplu;
      rhomin = cells[0]._rho - bondi_density(rmin_inv);
      rhoplu = cells[0]._rho - cells[1]._rho;
      umin = cells[0]._u - bondi_velocity(rmin_inv);
      uplu = cells[0]._u - cells[1]._u;
      if (istep < NSTEP_RELAX) {
        Pmin = cells[0]._P - bondi_pressure(rmin_inv);
      } else {
        Pmin = cells[0]._P - 100. * bondi_pressure(rmin_inv);
      }
      Pplu = cells[0]._P - cells[1]._P;
      if (std::abs(rhomin) < std::abs(rhoplu)) {
        cells[0]._grad_rho = rhomin / CELLSIZE;
      } else {
        cells[0]._grad_rho = rhoplu / CELLSIZE;
      }
      if (std::abs(umin) < std::abs(uplu)) {
        cells[0]._grad_u = umin / CELLSIZE;
      } else {
        cells[0]._grad_u = uplu / CELLSIZE;
      }
      if (std::abs(Pmin) < std::abs(Pplu)) {
        cells[0]._grad_P = Pmin / CELLSIZE;
      } else {
        cells[0]._grad_P = Pplu / CELLSIZE;
      }
    }

    // upper boundary
    {
      double rhomin, rhoplu, umin, uplu, Pmin, Pplu, dmin, dplu;
      rhomin = cells[NCELL + 1]._rho - cells[NCELL]._rho;
      rhoplu = cells[NCELL + 1]._rho - bondi_density(rmax_inv);
      umin = cells[NCELL + 1]._u - cells[NCELL]._u;
      uplu = cells[NCELL + 1]._u - bondi_velocity(rmax_inv);
      Pmin = cells[NCELL + 1]._P - cells[NCELL]._P;
      Pplu = cells[NCELL + 1]._P - bondi_pressure(rmax_inv);
      if (std::abs(rhomin) < std::abs(rhoplu)) {
        cells[NCELL + 1]._grad_rho = rhomin / CELLSIZE;
      } else {
        cells[NCELL + 1]._grad_rho = rhoplu / CELLSIZE;
      }
      if (std::abs(umin) < std::abs(uplu)) {
        cells[NCELL + 1]._grad_u = umin / CELLSIZE;
      } else {
        cells[NCELL + 1]._grad_u = uplu / CELLSIZE;
      }
      if (std::abs(Pmin) < std::abs(Pplu)) {
        cells[NCELL + 1]._grad_P = Pmin / CELLSIZE;
      } else {
        cells[NCELL + 1]._grad_P = Pplu / CELLSIZE;
      }
    }
#elif BOUNDARIES == BOUNDARIES_REFLECTIVE
    cells[0]._grad_rho = cells[1]._grad_rho;
    cells[0]._grad_u =
        -cells[1]._grad_u - cells[1]._grad_u -
        4. * cells[0]._u / (cells[0]._midpoint - cells[1]._midpoint);
    cells[0]._grad_P = cells[1]._grad_P;
    cells[NCELL + 1]._grad_rho = -cells[NCELL]._grad_rho;
    cells[NCELL + 1]._grad_u =
        -cells[NCELL]._grad_u -
        4. * cells[NCELL]._u /
            (cells[NCELL]._midpoint - cells[NCELL + 1]._midpoint);
    cells[NCELL + 1]._grad_P = -cells[NCELL]._grad_P;
#endif

#ifdef NO_GRADIENTS
    // reset all gradients to zero to disable the second order scheme
    for (unsigned int i = 0; i < NCELL + 2; ++i) {
      cells[i]._grad_rho = 0.;
      cells[i]._grad_u = 0.;
      cells[i]._grad_P = 0.;
    }
#endif

    // evolve all primitive variables forward in time for half a time step
    // using the Euler equations and the spatial gradients within the cells
    for (unsigned int i = 0; i < NCELL + 2; ++i) {
      const double rho = cells[i]._rho;
      const double u = cells[i]._u;
      const double P = cells[i]._P;
      cells[i]._rho -=
          0.5 * DT * (rho * cells[i]._grad_u + u * cells[i]._grad_rho);
      cells[i]._u -= 0.5 * DT * (u * cells[i]._grad_u + cells[i]._grad_P / rho);
      cells[i]._P -=
          0.5 * DT * (GAMMA * P * cells[i]._grad_u + u * cells[i]._grad_P);

#if POTENTIAL != POTENTIAL_NONE
      // add the gravitational time step prediction
      cells[i]._u += 0.5 * DT * cells[i]._a;
#endif
    }

    // do the flux exchange
    // we always do the flux for the right face of the cell, which means we
    // start from the left ghost and do not include the right ghost
    for (unsigned int i = 0; i < NCELL + 1; ++i) {
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

    // add the spherical source term (see above)
    for (unsigned int i = 1; i < NCELL + 1; ++i) {
      const double r = cells[i]._midpoint;
      const double Ui[3] = {cells[i]._m / cells[i]._V,
                            cells[i]._p / cells[i]._V,
                            cells[i]._E / cells[i]._V};
      double K1[3];
      K1[0] = -DT * Ui[1] / r;
      K1[1] = -DT * Ui[1] * Ui[1] / Ui[0] / r;
      const double p1 = (GAMMA - 1.) * (Ui[2] - 0.5 * Ui[1] * Ui[1] / Ui[0]);
      K1[2] = -DT * Ui[1] * (Ui[2] + p1) / Ui[0] / r;
      double U[3] = {Ui[0] + K1[0], Ui[1] + K1[1], Ui[2] + K1[2]};
      double K2[3];
      K2[0] = -DT * U[1] / r;
      K2[1] = -DT * U[1] * U[1] / U[0] / r;
      const double p2 = (GAMMA - 1.) * (U[2] - 0.5 * U[1] * U[1] / U[0]);
      K2[2] = -DT * U[1] * (U[2] + p2) / U[0] / r;
      U[0] = Ui[0] + 0.5 * (K1[0] + K2[0]);
      U[1] = Ui[1] + 0.5 * (K1[1] + K2[1]);
      U[2] = Ui[2] + 0.5 * (K1[2] + K2[2]);

      cells[i]._m = U[0] * cells[i]._V;
      cells[i]._p = U[1] * cells[i]._V;
      cells[i]._E = U[2] * cells[i]._V;
    }

// add gravitational acceleration (see above)
#if POTENTIAL == POTENTIAL_POINT_MASS
    for (unsigned int i = 1; i < NCELL + 1; ++i) {
      const double r = cells[i]._midpoint;
      const double a = -G * MASS_POINT_MASS / (r * r);
      cells[i]._a = a;
      const double m = cells[i]._V * cells[i]._rho;
      cells[i]._p += 0.5 * DT * a * m;
    }
#endif
  }

  // write the final snapshot
  write_snapshot(NSTEP, cells);
  write_binary_snapshot(cells);

  // all went well: return with exit code 0
  return 0;
}
