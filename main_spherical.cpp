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
 * @file main_spherical.cpp
 *
 * @brief Main program.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */

// project includes
#include "Bondi.hpp"             // for EOS_BONDI, BOUNDARIES_BONDI, IC_BONDI
#include "Boundaries.hpp"        // for non Bondi boundary conditions
#include "Cell.hpp"              // Cell class
#include "EOS.hpp"               // for non Bondi equations of state
#include "HLLCRiemannSolver.hpp" // fast HLLC Riemann solver
#include "IC.hpp"                // general initial condition interface
#include "LogFile.hpp"           // log file output
#include "Potential.hpp"         // external gravity
#include "RiemannSolver.hpp"     // slow exact Riemann solver
#include "SafeParameters.hpp"    // safe way to include Parameter.hpp
#include "Spherical.hpp"         // spherical source terms
#include "Timer.hpp"             // program timers
#include "Units.hpp"             // unit information

// standard libraries
#include <cfloat>
#include <cmath>
#include <cstdint>
#include <ctime>
#include <fstream>
#include <iostream>
#include <omp.h>
#include <sstream>

/*! @brief Activate this to disable fancy log output. */
#define NO_LOGFILE

/**
 * @brief Get the current time as a string.
 *
 * @return Current system time as a string with format YYYY:MM:DD:HH:MM:SS.
 */
static std::string get_timestamp() {
  const std::time_t timestamp = std::time(nullptr);
  const std::tm *time = std::localtime(&timestamp);
  std::stringstream timestream;
  timestream << (time->tm_year + 1900) << ":";
  if (time->tm_mon < 9) {
    timestream << "0";
  }
  timestream << (time->tm_mon + 1) << ":";
  if (time->tm_mday < 10) {
    timestream << "0";
  }
  timestream << time->tm_mday << ":";
  if (time->tm_hour < 10) {
    timestream << "0";
  }
  timestream << time->tm_hour << ":";
  if (time->tm_min < 10) {
    timestream << "0";
  }
  timestream << time->tm_min << ":";
  if (time->tm_sec < 10) {
    timestream << "0";
  }
  timestream << time->tm_sec;
  return timestream.str();
}

/**
 * @brief Write a snapshot with the given index.
 *
 * @param istep Index of the snapshot file.
 * @param cells Cells to write.
 * @param ncell Number of cells.
 */
void write_snapshot(uint_fast64_t istep, double time, const Cell *cells,
                    const unsigned int ncell) {
  std::stringstream filename;
  filename << "snapshot_";
  filename.fill('0');
  filename.width(4);
  filename << istep;
  filename << ".txt";
  std::cout << "Writing snapshot " << filename.str() << std::endl;
  std::ofstream ofile(filename.str().c_str());
  ofile << "# time: " << time * UNIT_TIME_IN_SI << "\n";
  for (uint_fast32_t i = 1; i < ncell + 1; ++i) {
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
 * @param cells Cells to write.
 * @param ncell Number of cells.
 */
void write_binary_snapshot(const Cell *cells, const unsigned int ncell) {
  std::ofstream ofile("lastsnap.dat");
  for (uint_fast32_t i = 1; i < ncell + 1; ++i) {
    ofile.write(reinterpret_cast<const char *>(&cells[i]._rho), sizeof(double));
    ofile.write(reinterpret_cast<const char *>(&cells[i]._u), sizeof(double));
    ofile.write(reinterpret_cast<const char *>(&cells[i]._P), sizeof(double));
    ofile.write(reinterpret_cast<const char *>(&cells[i]._a), sizeof(double));
  }
}

/**
 * @brief Practical names for log entry numbers.
 */
enum LogEntry {
  LOGENTRY_DENSITY = 0,
  LOGENTRY_VELOCITY,
  LOGENTRY_PRESSURE,
  LOGENTRY_NFRAC,
  NUMBER_OF_LOGENTRIES
};

/**
 * @brief Check if the value for the given log entry variable has changed
 * significantly since the last output.
 *
 * @param logentry Log entry identifier.
 * @param cell Cell to check.
 */
static inline bool changed(const int logentry, const Cell &cell) {
  // tolerance: If the relative difference of the value and the last outputted
  // value is less than this value, no output is written
  // Should probably become a parameter at some point...
  static const double tol = 1.e-3;
  switch (logentry) {
  case LOGENTRY_DENSITY:
    return std::abs(cell._rho - cell._last_rho) >
           tol * std::abs(cell._rho + cell._last_rho);
  case LOGENTRY_VELOCITY:
    return std::abs(cell._u - cell._last_u) >
           tol * std::abs(cell._u + cell._last_u);
  case LOGENTRY_PRESSURE:
    return std::abs(cell._P - cell._last_P) >
           tol * std::abs(cell._P + cell._last_P);
  case LOGENTRY_NFRAC:
    return std::abs(cell._nfac - cell._last_nfac) >
           tol * std::abs(cell._nfac + cell._last_nfac);
  default:
    return false;
  }
}

/**
 * @brief Get the value for the given log entry variable.
 *
 * @param logentry Log entry identifier.
 * @param cell Cell for which we want the variable.
 * @return Value of the variable.
 */
static inline double get_value(const int logentry, Cell &cell) {
  switch (logentry) {
  case LOGENTRY_DENSITY:
    cell._last_rho = cell._rho;
    return cell._rho * UNIT_DENSITY_IN_SI;
  case LOGENTRY_VELOCITY:
    cell._last_u = cell._u;
    return cell._u * UNIT_VELOCITY_IN_SI;
  case LOGENTRY_PRESSURE:
    cell._last_P = cell._P;
    return cell._P * UNIT_PRESSURE_IN_SI;
  case LOGENTRY_NFRAC:
    cell._last_nfac = cell._nfac;
    return cell._nfac;
  default:
    return 0.;
  }
}

/**
 * @brief Write significantly changed variables to the log file.
 *
 * @param log LogFile to write to.
 * @param cells Cells to write.
 * @param ncell Number of cells.
 * @param time Current simulation time (in internal units of T).
 * @param full_dump If set to True, dumps all cells irrespective of variable
 * changes.
 */
static inline void write_logfile(LogFile &log, Cell *cells,
                                 const unsigned int ncell, const double time,
                                 bool full_dump = false) {
#ifndef NO_LOGFILE
  if (full_dump) {
    // full dump: write all particles
    for (uint_fast16_t i = 1; i < ncell + 1; ++i) {
      for (int logentry = 0; logentry < NUMBER_OF_LOGENTRIES; ++logentry) {
        unsigned long previous_entry = cells[i]._last_entry;
        cells[i]._last_entry = log.get_current_position();
        previous_entry = cells[i]._last_entry - previous_entry;
        log.write(previous_entry);
        log.write(cells[i]._index);
        log.write(logentry);
        log.write(time);
        log.write(get_value(logentry, cells[i]));
      }
    }
  } else {
    // only write interesting quantities
    for (uint_fast16_t i = 1; i < ncell + 1; ++i) {
      for (int logentry = 0; logentry < NUMBER_OF_LOGENTRIES; ++logentry) {
        if (changed(logentry, cells[i])) {
          unsigned long previous_entry = cells[i]._last_entry;
          cells[i]._last_entry = log.get_current_position();
          previous_entry = cells[i]._last_entry - previous_entry;
          log.write(previous_entry);
          log.write(cells[i]._index);
          log.write(logentry);
          log.write(time);
          log.write(get_value(logentry, cells[i]));
        }
      }
    }
  }
#endif // NO_LOGFILE
}

/**
 * @brief Round the given integer down to the nearest power of 2.
 *
 * @param x Integer.
 * @return Nearest lower power of 2.
 */
static inline uint_fast64_t round_power2_down(uint_fast64_t x) {
  if (x == 0) {
    std::cerr << "Zero time step!" << std::endl;
    exit(1);
  }
  --x;
  x |= (x >> 1);
  x |= (x >> 2);
  x |= (x >> 4);
  x |= (x >> 8);
  x |= (x >> 16);
  x |= (x >> 32);
  x >>= 1;
  ++x;
  if (x == 0) {
    std::cerr << "Zero time step!" << std::endl;
    exit(1);
  }
  return x;
}

/**
 * @brief Get the conserved energy for a shell based on the corresponding 1D
 * cell properties.
 *
 * This function corrects for the difference between the 1D cell energy and the
 * 3D shell energy.
 *
 * @param cell Cell.
 * @return Shell energy.
 */
static inline double get_shell_energy(const Cell &cell) {
  const double rho = cell._rho;
  const double u = cell._u;
  const double P = cell._P;
#if DIMENSIONALITY == DIMENSIONALITY_1D
  const double Vs = cell._V;
#else
  const double Ru = cell._uplim;
  const double Rl = cell._lowlim;
  const double Vs = 4. * M_PI / 3. * (Ru * Ru * Ru - Rl * Rl * Rl);
#endif
  const double E = P * Vs / (GAMMA - 1.) + 0.5 * rho * Vs * u * u;
  return E;
}

/**
 * @brief Main simulation program.
 *
 * Usage: ./HydroCodeSpherical1D [ncell [ic_file_name [transition_width
 *        [bondi_pressure_contrast] ] ] ]
 * Valid command line arguments are:
 *  - ncell: Number of cells to use.
 *  - ic_file_name: Name of the initial condition file (if IC_FILE was selected
 *    when configuring the code). Note that the number of values in the file
 *    should match the value of "ncell".
 *  - transition_width: Width of the linear transition region between ionised
 *    and neutral region (if IONISATION_TRANSITION_SMOOTH was selected when
 *    configuring the code). Should be given in internal length units.
 *  - bondi_pressure_contrast: Pressure contrast between ionised and neutral
 *    region. The pressure contrast consist of the ratio of ionised and neutral
 *    temperature AND the factor two change in mean particle mass between
 *    ionised and neutral gas.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  // time the program
  Timer total_time;
  total_time.start();

  // initialize the optional parameters with their default values
  // default values are given in Parameters.hpp.in, DerivedParameters.hpp and
  // Bondi.hpp
  unsigned int ncell = NCELL;
  std::string ic_file_name(IC_FILE_NAME);
  double transition_width = IONISATION_TRANSITION_WIDTH;
  double bondi_pressure_contrast = BONDI_PRESSURE_CONTRAST;
  // disable unused variable warnings
  (void)bondi_pressure_contrast;

  // now overwrite with the actual command line parameters (if specified)
  if (argc > 1) {
    ncell = atoi(argv[1]);
  }
  if (argc > 2) {
    ic_file_name = argv[2];
  }
  if (argc > 3) {
    transition_width = atof(argv[3]);
  }
  if (argc > 4) {
    bondi_pressure_contrast = atof(argv[4]);
  }

  // output: most of this was useful at some point
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
            << G_INTERNAL *
                   (UNIT_LENGTH_IN_SI * UNIT_LENGTH_IN_SI * UNIT_LENGTH_IN_SI /
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
  std::cout << "Neutral Bondi radius: " << RBONDI << " ("
            << RBONDI * UNIT_LENGTH_IN_SI / AU_IN_SI << " AU)" << std::endl;
  std::cout << "Density at R_Bondi: "
            << bondi_density(RBONDI * UNIT_LENGTH_IN_SI / (20. * AU_IN_SI)) *
                   UNIT_DENSITY_IN_SI
            << std::endl;
#endif

  std::cout << "Initial ionisation radius: "
            << INITIAL_IONISATION_RADIUS * UNIT_LENGTH_IN_SI / AU_IN_SI
            << " AU (" << INITIAL_IONISATION_RADIUS << ")" << std::endl;

  std::cout << "Point mass: " << MASS_POINT_MASS * UNIT_MASS_IN_SI << " kg"
            << std::endl;

  std::cout << "Useful units:" << std::endl;
  std::cout << "Point mass: " << MASS_POINT_MASS * UNIT_MASS_IN_MSOL << " Msol"
            << std::endl;
  std::cout << "Total simulation time: " << MAXTIME * UNIT_TIME_IN_YR << " yr"
            << std::endl;
  std::cout << "Time in between snapshots: "
            << (MAXTIME / NUMBER_OF_SNAPS) * UNIT_TIME_IN_YR << " yr "
            << std::endl;
  std::cout << "Minimum radius: " << RMIN * UNIT_LENGTH_IN_AU << " AU (" << RMIN
            << ")" << std::endl;
  std::cout << "Maximum radius: " << RMAX * UNIT_LENGTH_IN_AU << " AU (" << RMAX
            << ")" << std::endl;

// figure out how many threads we are using and tell the user about this
#pragma omp parallel
  {
#pragma omp single
    {
      int num_thread = omp_get_num_threads();
      std::cout << "Running on " << num_thread << " thread(s)." << std::endl;
    }
  }

  // initialize the time line used for time stepping
  // we use a classical power of 2 integer time line as e.g. Gadget2
  const double maxtime = MAXTIME;
  const uint_fast64_t integer_maxtime = 0x8000000000000000; // 2^63
  const double time_conversion_factor = maxtime / integer_maxtime;

  // create the 1D spherical grid
  // we create 2 ghost cells to the left and to the right of the simulation box
  // to handle boundary conditions
  Cell *cells = new Cell[ncell + 2];
#pragma omp parallel for
  for (uint_fast32_t i = 0; i < ncell + 2; ++i) {
    // cell positions (lower limit, center and upper limit) are precomputed for
    // maximal efficiency
    cells[i]._lowlim = RMIN + (i - 1.) * CELLSIZE;
    cells[i]._midpoint = RMIN + (i - 0.5) * CELLSIZE;
    cells[i]._uplim = RMIN + i * CELLSIZE;
    cells[i]._V = CELLSIZE;
    cells[i]._integer_dt = 0;
    // initialize the time step to a sensible value: the requested snapshot time
    // interval
    cells[i]._dt = (MAXTIME / NUMBER_OF_SNAPS);
    // only actual cells have an index in the range [0, ncell[. The two ghost
    // cells are excluded by the bit of conditional magic below.
    cells[i]._index = (i != 0 && i != ncell + 2) ? (i - 1) : ncell + 2;
  }

  // set up the initial condition
  // this bit is handled by IC.hpp, and specific implementations in ICFile.hpp
  // (if configured with IC_FILE), Bondi.hpp (if configured with IC_BONDI), or
  // Sod.hpp (if configured with IC_SOD).
  initialize(cells, ncell);

  // Courant factor for the CFL time step criterion
  // we use a very conservative value
  const double courant_factor = COURANT_FACTOR;

  std::cout << "Courant factor: " << courant_factor << std::endl;

  // initialize snapshot variables
  const uint_fast64_t snaptime = integer_maxtime / NUMBER_OF_SNAPS;
  uint_fast64_t isnap = 0;

  // convert the input primitive variables into conserved variables, and compute
  // the initial time step
  // we use a global time step, which is the minimum time step among all cells
  uint_fast64_t min_integer_dt = snaptime;
  double Etot = 0.;
#pragma omp parallel for reduction(min : min_integer_dt) reduction(+ : Etot)
  // convert primitive variables to conserved variables
  for (uint_fast32_t i = 1; i < ncell + 1; ++i) {
    // apply the equation of state to get the initial pressure (if necessary)
    // this bit is handled by EOS.hpp and Bondi.hpp (for EOS_BONDI)
    initial_pressure(cells[i]);

    // use the cell volume to convert primitive into conserved variables
    cells[i]._m = cells[i]._rho * cells[i]._V;
    cells[i]._p = cells[i]._m * cells[i]._u;
    cells[i]._E = cells[i]._P * cells[i]._V / (GAMMA - 1.) +
                  0.5 * cells[i]._u * cells[i]._p;

    Etot += get_shell_energy(cells[i]);

    // time step criterion
    const double cs =
        std::sqrt(GAMMA * cells[i]._P / cells[i]._rho) + std::abs(cells[i]._u);
    const double dt = courant_factor * cells[i]._V / cs;
    const uint_fast64_t integer_dt =
        (dt < maxtime) ? (dt / maxtime) * integer_maxtime : integer_maxtime;
    min_integer_dt = std::min(min_integer_dt, integer_dt);
    if (min_integer_dt == 0) {
      std::cerr << "Cell pushes time step to 0!" << std::endl;
      std::cerr << "dt: " << dt << std::endl;
      std::cerr << "maxtime: " << maxtime << std::endl;
      std::cerr << "integer_dt: " << integer_dt << std::endl;
      std::cerr << "explicitly: " << (dt / maxtime) * integer_maxtime
                << std::endl;
      std::cerr << "cs: " << cs << std::endl;
      cells[i].print();
      exit(1);
    }

    // initialize variables used for the log file
    cells[i]._last_rho = cells[i]._rho;
    cells[i]._last_u = cells[i]._u;
    cells[i]._last_P = cells[i]._P;
    cells[i]._last_nfac = cells[i]._nfac;
    cells[i]._last_entry = 0;
  }

  std::cout << "Initial minimum timestep: " << min_integer_dt << std::endl;
  std::cout << "Initial total energy: " << Etot * UNIT_ENERGY_IN_SI << " J"
            << std::endl;

  std::ofstream energy_file("energy.txt");
  energy_file << "0.\t" << Etot * UNIT_ENERGY_IN_SI << "\n";

  // initialize the log file and write the first entry
  LogFile logfile("logfile.dat", 100);
  write_logfile(logfile, cells, ncell, 0., true);

  // set cell time steps
  // round min_integer_dt to closest smaller power of 2
  uint_fast64_t global_integer_dt = round_power2_down(min_integer_dt);
  std::cout << "Rounded minimum timestep: " << global_integer_dt << std::endl;
#pragma omp parallel for
  for (uint_fast32_t i = 0; i < ncell + 2; ++i) {
    cells[i]._integer_dt = global_integer_dt;
    cells[i]._dt = cells[i]._integer_dt * time_conversion_factor;
  }

  // initialize boundary condition and ionisation variables
  // these bits are handled in EOS.hpp (and Bondi.hpp for EOS_BONDI), and
  // Boundaries.hpp (and Bondi.hpp for BOUNDARIES_BONDI).
  boundary_conditions_initialize();
  ionisation_initialize();

// initialize the Riemann solver
// we use a fast HLLC solver
// replace "HLLCRiemannSolver" with "RiemannSolver" to use a slower, exact
// solver
#if RIEMANNSOLVER_TYPE == RIEMANNSOLVER_TYPE_HLLC
  HLLCRiemannSolver solver(GAMMA);
#elif RIEMANNSOLVER_TYPE == RIEMANNSOLVER_TYPE_EXACT
  RiemannSolver solver(GAMMA);
#endif

  // initialize some variables used to guesstimate the remaing run time
  Timer progress_timer;
  Timer step_time;
  double time_since_last = 0.;
  double time_since_start = 0.;
  unsigned int steps_since_last = 0;
  // initialize the time stepping
  uint_fast64_t current_integer_time = 0;
  uint_fast64_t current_integer_dt = global_integer_dt;
  std::cout << "Initial system time step: "
            << current_integer_dt * time_conversion_factor << std::endl;
  std::cout << "Maximum time: " << integer_maxtime * time_conversion_factor
            << std::endl;
  // main simulation loop: perform NSTEP steps
  while (current_integer_time < integer_maxtime) {

    // start the step timer
    step_time.start();

    // add the spherical source term. Handled by Spherical.hpp
    add_spherical_source_term();

    // do first gravity kick, handled by Potential.hpp
    do_gravity();

    // do ionisation, handled by EOS.hpp (and Bondi.hpp for EOS_BONDI).
    do_ionisation();

    // update the primitive variables based on the values of the conserved
    // variables and the current cell volume
    // also compute the new time step
    min_integer_dt = snaptime;
    Etot = 0.;
#pragma omp parallel for reduction(min : min_integer_dt) reduction(+ : Etot)
    for (uint_fast32_t i = 1; i < ncell + 1; ++i) {
      cells[i]._rho = cells[i]._m / cells[i]._V;
      cells[i]._u = cells[i]._p / cells[i]._m;
      // the pressure update depends on the equation of state
      // this is handled in EOS.hpp (and Bondi.hpp for EOS_BONDI)
      update_pressure(cells[i]);

      Etot += get_shell_energy(cells[i]);

      const double cs = std::sqrt(GAMMA * cells[i]._P / cells[i]._rho) +
                        std::abs(cells[i]._u);
      const double dt = courant_factor * cells[i]._V / cs;
      const uint_fast64_t integer_dt =
          (dt < maxtime) ? (dt / maxtime) * integer_maxtime : integer_maxtime;
      min_integer_dt = std::min(min_integer_dt, integer_dt);
    }

    energy_file << current_integer_time * time_conversion_factor << "\t"
                << Etot * UNIT_ENERGY_IN_SI << "\n";

    // now is the time to write to the log file
    write_logfile(logfile, cells, ncell,
                  current_integer_time * time_conversion_factor);

    // round min_integer_dt to closest smaller power of 2
    global_integer_dt = round_power2_down(min_integer_dt);
    // make sure the remaining time can be filled *exactly* with the current
    // time step
    while ((integer_maxtime - current_integer_time) % global_integer_dt > 0) {
      global_integer_dt >>= 1;
    }
// set the new cell time steps
#pragma omp parallel for
    for (uint_fast32_t i = 0; i < ncell + 2; ++i) {
      cells[i]._integer_dt = global_integer_dt;
      cells[i]._dt = cells[i]._integer_dt * time_conversion_factor;
    }
    current_integer_dt = global_integer_dt;

    // check if it is time for a status update
    if (progress_timer.interval() >= STATUS_UPDATE_INTERVAL) {
      // yes: display some statistics and a guesstimate of the remaining run
      // time
      const double pct = current_integer_time * 100. / integer_maxtime;
      std::cout << get_timestamp() << "\t" << ncell << ": "
                << "time " << current_integer_time * time_conversion_factor
                << " of " << maxtime << " (" << pct << " %)" << std::endl;
      std::cout << "\t\t\tSystem time step: "
                << current_integer_dt * time_conversion_factor << std::endl;
      const double avg_time_since_last = time_since_last / steps_since_last;
      std::cout << "\t\t\tAverage time per step: " << avg_time_since_last
                << " s" << std::endl;
      time_since_start += time_since_last;
      const double time_to_go = time_since_start * (100. - pct) / pct;
      std::cout << "\t\t\tEstimated time to go: " << time_to_go << " s"
                << std::endl;
#if EOS == EOS_BONDI
      // we added this bit for the case where we want to add accreted material
      // to the central mass (currently not used)
      std::cout << "\t\t\tCentral mass: " << central_mass << " ("
                << (central_mass / MASS_POINT_MASS) << ")" << std::endl;
#endif
      std::cout << "Total energy: " << Etot * UNIT_ENERGY_IN_SI << " J"
                << std::endl;
      // reset guesstimate counters
      time_since_last = 0.;
      steps_since_last = 0;
      // reset the progress timer
      progress_timer.start();
    }

    // check if we need to output a snapshot
    if (current_integer_time >= isnap * snaptime) {
      // write the actual snapshot
      write_snapshot(isnap, current_integer_time * time_conversion_factor,
                     cells, ncell);
      ++isnap;
    }

    // apply boundary conditions
    // handled by Boundaries.hpp (and Bondi.hpp for BOUNDARIES_BONDI)
    boundary_conditions_primitive_variables();

// compute slope limited gradients for the primitive variables in each cell
#pragma omp parallel for
    for (uint_fast32_t i = 1; i < ncell + 1; ++i) {
      const double dx = cells[i + 1]._midpoint - cells[i - 1]._midpoint;
      const double dx_inv = 1. / dx;
      const double half_dx = 0.5 * dx;

      const double gradrho = (cells[i + 1]._rho - cells[i - 1]._rho) * dx_inv;
      const double rhomax = std::max(cells[i - 1]._rho, cells[i + 1]._rho);
      const double rhomin = std::min(cells[i - 1]._rho, cells[i + 1]._rho);
      const double rho_ext_plu = half_dx * gradrho;
      const double rho_ext_min = -half_dx * gradrho;
      const double rhoextmax = std::max(rho_ext_min, rho_ext_plu);
      const double rhoextmin = std::min(rho_ext_min, rho_ext_plu);
      const double alpha_rho =
          (gradrho != 0.)
              ? std::min(1.,
                         0.5 * std::min((rhomax - cells[i]._rho) / rhoextmax,
                                        (rhomin - cells[i]._rho) / rhoextmin))
              : 1.;
      cells[i]._grad_rho = alpha_rho * gradrho;

      const double gradu = (cells[i + 1]._u - cells[i - 1]._u) * dx_inv;
      const double umax = std::max(cells[i - 1]._u, cells[i + 1]._u);
      const double umin = std::min(cells[i - 1]._u, cells[i + 1]._u);
      const double u_ext_plu = half_dx * gradu;
      const double u_ext_min = -half_dx * gradu;
      const double uextmax = std::max(u_ext_min, u_ext_plu);
      const double uextmin = std::min(u_ext_min, u_ext_plu);
      const double alpha_u =
          (gradu != 0.)
              ? std::min(1., 0.5 * std::min((umax - cells[i]._u) / uextmax,
                                            (umin - cells[i]._u) / uextmin))
              : 1.;
      cells[i]._grad_u = alpha_u * gradu;

      const double gradP = (cells[i + 1]._P - cells[i - 1]._P) * dx_inv;
      const double Pmax = std::max(cells[i - 1]._P, cells[i + 1]._P);
      const double Pmin = std::min(cells[i - 1]._P, cells[i + 1]._P);
      const double P_ext_plu = half_dx * gradP;
      const double P_ext_min = -half_dx * gradP;
      const double Pextmax = std::max(P_ext_min, P_ext_plu);
      const double Pextmin = std::min(P_ext_min, P_ext_plu);
      const double alpha_P =
          (gradP != 0.)
              ? std::min(1., 0.5 * std::min((Pmax - cells[i]._P) / Pextmax,
                                            (Pmin - cells[i]._P) / Pextmin))
              : 1.;
      cells[i]._grad_P = alpha_P * gradP;
    }

    // apply boundary conditions for the gradients
    // handled by Boundaries.hpp (and Bondi.hpp for BOUNDARIES_BONDI)
    boundary_conditions_gradients();

#if HYDRO_ORDER == 1
// reset all gradients to zero to disable the second order scheme
#pragma omp parallel for
    for (uint_fast32_t i = 0; i < ncell + 2; ++i) {
      cells[i]._grad_rho = 0.;
      cells[i]._grad_u = 0.;
      cells[i]._grad_P = 0.;
    }
#endif

// evolve all primitive variables forward in time for half a time step
// using the Euler equations and the spatial gradients within the cells
#pragma omp parallel for
    for (uint_fast32_t i = 0; i < ncell + 2; ++i) {
      const double half_dt = 0.5 * cells[i]._dt;
      const double rho = cells[i]._rho;
      const double u = cells[i]._u;
      const double P = cells[i]._P;
      cells[i]._rho -=
          half_dt * (rho * cells[i]._grad_u + u * cells[i]._grad_rho);
      if (rho > 0.) {
        cells[i]._u -=
            half_dt * (u * cells[i]._grad_u + cells[i]._grad_P / rho);
      }
      cells[i]._P -=
          half_dt * (GAMMA * P * cells[i]._grad_u + u * cells[i]._grad_P);

      if (cells[i]._rho < 0.) {
        cells[i]._rho = rho;
      }
      if (cells[i]._P < 0.) {
        cells[i]._P = P;
      }

      // add gravity prediction. Handled by Potential.hpp.
      add_gravitational_prediction(cells[i], half_dt);

      cells[i]._left_flux[0] = 0.;
      cells[i]._left_flux[1] = 0.;
      cells[i]._left_flux[2] = 0.;
      cells[i]._right_flux[0] = 0.;
      cells[i]._right_flux[1] = 0.;
      cells[i]._right_flux[2] = 0.;
    }

    // compute the fluxes
    // to avoid thread concurrency, we first compute the fluxes in separate
    // variables and then update the conserved quantities in a second loop
#pragma omp parallel for
    for (uint_fast32_t i = 1; i < ncell + 2; ++i) {
      const double dt = cells[i]._dt;
      // get the variables in the left and right state
      const double rhoL = cells[i - 1]._rho;
      const double uL = cells[i - 1]._u;
      const double PL = cells[i - 1]._P;
      const double rhoR = cells[i]._rho;
      const double uR = cells[i]._u;
      const double PR = cells[i]._P;

      // do the second order spatial reconstruction
      const double dmin = 0.5 * (cells[i]._midpoint - cells[i - 1]._midpoint);
      const double dplu = -dmin;
      double rhoL_dash, rhoR_dash, uL_dash, uR_dash, PL_dash, PR_dash;
      rhoL_dash = rhoL + dmin * cells[i - 1]._grad_rho;
      uL_dash = uL + dmin * cells[i - 1]._grad_u;
      PL_dash = PL + dmin * cells[i - 1]._grad_P;
      rhoR_dash = rhoR + dplu * cells[i]._grad_rho;
      uR_dash = uR + dplu * cells[i]._grad_u;
      PR_dash = PR + dplu * cells[i]._grad_P;

      if (i == 1) {
        boundary_left(rhoL_dash, uL_dash, PL_dash, rhoR_dash, uR_dash, PR_dash);
      }
      if (i == ncell) {
        boundary_right(rhoL_dash, uL_dash, PL_dash, rhoR_dash, uR_dash,
                       PR_dash);
      }

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

      // solve the Riemann problem at the interface between the two cells
      double mflux, pflux, Eflux;
      solver.solve_for_flux(rhoL_dash, uL_dash, PL_dash, rhoR_dash, uR_dash,
                            PR_dash, mflux, pflux, Eflux);

      // set the left and right fluxes
      // (unless the corresponding cell is a ghost)
      if (i < ncell + 1) {
        cells[i]._left_flux[0] = dt * mflux;
        cells[i]._left_flux[1] = dt * pflux;
        cells[i]._left_flux[2] = dt * Eflux;
      } else {
        cells[i]._left_flux[0] = 0.;
        cells[i]._left_flux[1] = 0.;
        cells[i]._left_flux[2] = 0.;
      }
      if (i > 1) {
        cells[i - 1]._right_flux[0] = -dt * mflux;
        cells[i - 1]._right_flux[1] = -dt * pflux;
        cells[i - 1]._right_flux[2] = -dt * Eflux;
      } else {
        cells[i - 1]._right_flux[0] = 0.;
        cells[i - 1]._right_flux[1] = 0.;
        cells[i - 1]._right_flux[2] = 0.;
      }

      // call a special function for flux that crosses the inner outflow
      // boundary. This currently does not do anything.
      if (i == 1) {
        flux_into_inner_mask(dt * mflux);
      }
    }

    // now actually do the flux exchanges by updating the conserved variables
#pragma omp parallel for
    for (uint_fast32_t i = 1; i < ncell + 1; ++i) {
      cells[i]._m += cells[i]._left_flux[0] + cells[i]._right_flux[0];
      cells[i]._p += cells[i]._left_flux[1] + cells[i]._right_flux[1];
      cells[i]._E += cells[i]._left_flux[2] + cells[i]._right_flux[2];
    }

    // add the spherical source term
    // handled by Spherical.hpp
    add_spherical_source_term();

    // do the second gravity kick
    // handled by Potential.hpp
    do_gravity();

    // stop the step timer, and update guesstimate counters
    step_time.stop();
    time_since_last += step_time.value();
    ++steps_since_last;
    step_time.reset();

    // update the system time
    current_integer_time += current_integer_dt;
  }
  std::cout << "Final total energy: " << Etot * UNIT_ENERGY_IN_SI << " J"
            << std::endl;

  // write the final logfile entry
  write_logfile(logfile, cells, ncell,
                current_integer_time * time_conversion_factor, true);
  // close the log file
  logfile.close_file();

  // write the final snapshots
  write_snapshot(isnap, current_integer_time * time_conversion_factor, cells,
                 ncell);
  write_binary_snapshot(cells, ncell);

  // clean up: free cell memory
  delete[] cells;

  // stop timing the program and display run time information
  total_time.stop();
  std::cout << "Total program time: " << total_time.value() << " s."
            << std::endl;

  // all went well: return with exit code 0
  return 0;
}
