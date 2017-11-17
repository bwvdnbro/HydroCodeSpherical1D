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
#include "HLLCRiemannSolver.hpp"
#include "IC.hpp"
#include "Potential.hpp"
#include "RiemannSolver.hpp"
#include "SafeParameters.hpp"
#include "Spherical.hpp"
#include "Timer.hpp"
#include "Units.hpp"

// standard libraries
#include <cfloat>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <ctime>
#include <fcntl.h>
#include <fstream>
#include <iostream>
#include <omp.h>
#include <sstream>
#include <sys/mman.h>
#include <unistd.h>

//#define NO_LOGFILE

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
 * @param Cells to write.
 */
void write_binary_snapshot(const Cell *cells, const unsigned int ncell) {
  std::ofstream ofile("lastsnap_sod.dat");
  for (uint_fast32_t i = 1; i < ncell + 1; ++i) {
    ofile.write(reinterpret_cast<const char *>(&cells[i]._rho), sizeof(double));
    ofile.write(reinterpret_cast<const char *>(&cells[i]._u), sizeof(double));
    ofile.write(reinterpret_cast<const char *>(&cells[i]._P), sizeof(double));
    ofile.write(reinterpret_cast<const char *>(&cells[i]._a), sizeof(double));
  }
}

/**
 * @brief Custom version of the SWIFT dump struct.
 */
class LogFile {
public:
  /* The memory-mapped data of this dump. */
  void *_data;

  /* The size of the memory-mapped data, in bytes. */
  size_t _size;

  /* The number of bytes that have been dumped. */
  size_t _count;

  /* The offset of the data within the current file. */
  size_t _file_offset;

  /* The file with which this memory is associated. */
  int _fd;

  /* Mask containing the significant bits for page addresses. */
  size_t _page_mask;
};

/**
 * @brief Convert the given size in MB to a size in bytes.
 *
 * @param size_in_MB Size in MB.
 * @return Size in bytes.
 */
static inline size_t MB_to_bytes(size_t size_in_MB) { return size_in_MB << 20; }

/**
 * @brief Initialize the given LogFile instance.
 *
 * @param log LogFile to initialize.
 * @param filename Name of the log file in the file system.
 * @param size Size of the log file in memory. The log file stored on disk can
 * be (significantly) larger, but the part that is directly stored in memory
 * will be this size. Note that the actually allocated memory size might be
 * larger than this, since we can only memory-map regions with a system
 * dependent size.
 */
static inline void initialize_log(LogFile &log, std::string filename,
                                  size_t size) {
  // create the file
  // O_CREAT: create the file if it does not exist
  // O_RDWR: read/write access
  // S_IRUSR/IWUSR/IRGRP/IWGRP: set read/write access for both the user and
  //  group when creating the file (code 0660 in the SWIFT code)
  log._fd = open(filename.c_str(), O_CREAT | O_RDWR,
                 S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP);
  // get the page size and convert it into a page mask
  // e.g. on a system with page size 4096:
  //  sysconf(_SC_PAGE_SIZE) = ..0001 0000 0000 0000
  //  ANS - 1                = ..0000 1111 1111 1111
  //  ~ANS = page_mask       = ..1111 0000 0000 0000
  const size_t page_mask = ~(sysconf(_SC_PAGE_SIZE) - 1);
  // make sure the requested size is a multiple of the page size (as this is
  // required for mmap()
  // e.g. size = 5
  //  size + ~page_mask = 101 + 1111 1111 1111
  //                    = 1 0000 0000 0100
  //  (ANS) & mask      = 1 0000 0000 0000
  //                    = page_size
  // e.g. size = page_size + 2
  //  size + ~page_mask = 1 0000 0000 0010 + 1111 1111 1111
  //                    = 10 0000 0000 0010
  //  (ANS) & mask      = 10 0000 0000 0000
  //                    = 2 * page_size
  size = (size + ~page_mask) & page_mask;
  // make sure we have enough space on disk to store the memory mapped file
  // this preallocates the file size without actually writing anything
  // we tell the call we want all file offsets between 0 and size to be
  // available after this call
  posix_fallocate(log._fd, 0, size);
  // memory-map the file: make a synced copy of some part of the file in active
  // memory that we can directly access as if it was normal memory
  // parameters:
  //  nullptr: let the system decide where to allocate the memory
  //  size: size of the buffer we want to memory-map (should be a multiple of
  //   page_size!)
  //  PROT_WRITE: we want to write to the buffer
  //  MAP_SHARED: we want to share the memory with all other processes: every
  //   write to this region is synced in the file and directly reflected in all
  //   other processes that memory map the same region
  //  log._fd: this is the file we want to memory-map
  //  0: we want to start memory mapping from offset 0 in the file
  log._data = mmap(nullptr, size, PROT_WRITE, MAP_SHARED, log._fd, 0);
  // initialize the size of the memory mapped region, the number of bytes
  // already written to it, the current offset in the file and the page_mask
  // (which is used later to increase the file size)
  log._size = size;
  log._count = 0;
  log._file_offset = 0;
  log._page_mask = page_mask;
}

/**
 * @brief Increase the size of the given log file.
 *
 * This will grow the size of the file on disk and will reset the part of the
 * file that is directly mapped in memory.
 *
 * @param log LogFile.
 */
static inline void resize_log(LogFile &log) {
  // truncate the number of bytes already written to the file to a multiple of
  // the page size
  // we will shift the memory-mapped region to this point
  const size_t trunc_count = log._count & log._page_mask;
  // unmap the part of the memory-mapped region before trunc_count
  // if nothing was written to the file, we provide 1, which is rounded up to
  // a single page size
  munmap(log._data, trunc_count > 0 ? trunc_count : 1);
  // increase the file offset, which indicates the offset of the memory-mapped
  // region within the file on disk
  log._file_offset += trunc_count;
  // subtract the unmapped region from the counter that counts the bytes already
  // written to the mapped region
  log._count -= trunc_count;
  // make sure we have enough disk space for the next bit of log file
  posix_fallocate(log._fd, log._file_offset, log._size);
  // map the new bit of file
  // parameters are the same as above, but now we start from a later offset in
  // the file
  log._data = mmap(nullptr, log._size, PROT_WRITE, MAP_SHARED, log._fd,
                   log._file_offset);
}

/**
 * @brief Force synchronization between the memory-mapped log file and the
 * contents of the log file on disk.
 *
 * @param log LogFile.
 */
static inline void force_log_sync(LogFile &log) {
  // make sure the memory-mapped region is synced on disk
  // parameters:
  //  log._data: pointer to the beginning of the region we want synced
  //  log._count: number of bytes after that pointer we want to be synced
  //  MS_SYNC: force actual writing to disk and wait until it is finished
  msync(log._data, log._count, MS_SYNC);
}

/**
 * @brief Close the log file.
 *
 * This method also ensures the file has its definitive file size on disk (more
 * space might be allocated then necessary during the run).
 *
 * @param log LogFile.
 */
static inline void close_log(LogFile &log) {
  // unmap the actively memory-mapped region of the file
  munmap(log._data, log._count);
  // shrink the file to its actual size
  ftruncate(log._fd, log._file_offset + log._count);
  // close the file
  close(log._fd);
}

/**
 * @brief Write the given value to the given log file.
 *
 * @param log LogFile.
 * @param value Value to write.
 */
template <typename _type_>
static inline void write_log(LogFile &log, _type_ value) {
  const size_t vsize = sizeof(_type_);
  const size_t offset = log._count;
  log._count += vsize;
  memcpy(log._data + offset, &value, vsize);
}

enum LogEntry {
  LOGENTRY_DENSITY = 0,
  LOGENTRY_VELOCITY,
  LOGENTRY_PRESSURE,
  LOGENTRY_NFRAC,
  NUMBER_OF_LOGENTRIES
};

static inline bool changed(const int logentry, const Cell &cell) {
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

static inline void write_logfile(LogFile &log, Cell *cells,
                                 const unsigned int ncell, const double time,
                                 bool full_dump = false) {
#ifndef NO_LOGFILE
  if (full_dump) {
    // initial write: write all particles
    for (uint_fast16_t i = 1; i < ncell + 1; ++i) {
      for (int logentry = 0; logentry < NUMBER_OF_LOGENTRIES; ++logentry) {
        unsigned long previous_entry = cells[i]._last_entry;
        cells[i]._last_entry = log._count;
        previous_entry = cells[i]._last_entry - previous_entry;
        double timecopy = time;
        double value = get_value(logentry, cells[i]);
        write_log(log, previous_entry);
        write_log(log, cells[i]._index);
        write_log(log, logentry);
        write_log(log, timecopy);
        write_log(log, value);
      }
    }
  } else {
    // only write interesting quantities
    for (uint_fast16_t i = 1; i < ncell + 1; ++i) {
      for (int logentry = 0; logentry < NUMBER_OF_LOGENTRIES; ++logentry) {
        if (changed(logentry, cells[i])) {
          unsigned long previous_entry = cells[i]._last_entry;
          cells[i]._last_entry = log._count;
          previous_entry = cells[i]._last_entry - previous_entry;
          double timecopy = time;
          double value = get_value(logentry, cells[i]);
          write_log(log, previous_entry);
          write_log(log, cells[i]._index);
          write_log(log, logentry);
          write_log(log, timecopy);
          write_log(log, value);
        }
      }
    }
  }
#endif // NO_LOGFILE
}

static inline uint_fast64_t round_power2_down(uint_fast64_t x) {
  --x;
  x |= (x >> 1);
  x |= (x >> 2);
  x |= (x >> 4);
  x |= (x >> 8);
  x |= (x >> 16);
  x |= (x >> 32);
  x >>= 1;
  ++x;
  return x;
}

/**
 * @brief Main simulation program.
 *
 * @param argc Number of command line arguments (ignored).
 * @param argv Command line arguments (ignored).
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  Timer total_time;
  total_time.start();

  unsigned int ncell = NCELL;
  std::string ic_file_name(IC_FILE_NAME);
  double transition_width = IONIZATION_TRANSITION_WIDTH;
  double bondi_pressure_contrast = BONDI_PRESSURE_CONTRAST;

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
  std::cout << "Bondi radius: " << RBONDI << " ("
            << RBONDI * UNIT_LENGTH_IN_SI / AU_IN_SI << " AU)" << std::endl;
  std::cout << "Ionized Bondi radius: " << RBONDI_ION << " ("
            << RBONDI_ION * UNIT_LENGTH_IN_SI / AU_IN_SI << " AU)" << std::endl;
  std::cout << "Density at R_Bondi: "
            << bondi_density(RBONDI * UNIT_LENGTH_IN_SI / (20. * AU_IN_SI)) *
                   UNIT_DENSITY_IN_SI
            << std::endl;
#endif

  //  std::cout << "Ionizing luminosity: "
  //            << (4. * M_PI * 4.e-7 * bondi_Q(INITIAL_IONIZATION_RADIUS) *
  //                UNIT_DENSITY_IN_SI * UNIT_MASS_IN_SI / HYDROGEN_MASS_IN_SI /
  //                HYDROGEN_MASS_IN_SI)
  //            << std::endl;

  std::cout << "Initial ionization radius: "
            << INITIAL_IONIZATION_RADIUS * UNIT_LENGTH_IN_SI / AU_IN_SI
            << std::endl;

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

#pragma omp parallel
  {
#pragma omp single
    {
      int num_thread = omp_get_num_threads();
      std::cout << "Running on " << num_thread << " thread(s)." << std::endl;
    }
  }

  const double maxtime = DT * NSTEP;
  const uint_fast64_t integer_maxtime = 0x8000000000000000; // 2^63
  const double time_conversion_factor = maxtime / integer_maxtime;

  // create the 1D spherical grid
  // we create 2 ghost cells to the left and to the right of the simulation box
  // to handle boundary conditions
  Cell *cells = new Cell[ncell + 2];
#pragma omp parallel for
  for (uint_fast32_t i = 0; i < ncell + 2; ++i) {
    cells[i]._lowlim = RMIN + (i - 1.) * CELLSIZE;
    cells[i]._midpoint = RMIN + (i - 0.5) * CELLSIZE;
    cells[i]._uplim = RMIN + i * CELLSIZE;
    cells[i]._V = CELLSIZE;
    cells[i]._integer_dt = 0;
    cells[i]._dt = DT;
    cells[i]._index = (i != 0 && i != ncell + 2) ? (i - 1) : ncell + 2;
  }

  // set up the initial condition
  initialize(cells, ncell);

  const double courant_factor = 0.01;

  uint_fast64_t min_integer_dt = integer_maxtime;
#pragma omp parallel for reduction(min : min_integer_dt)
  // convert primitive variables to conserved variables
  for (uint_fast32_t i = 1; i < ncell + 1; ++i) {
    // apply the equation of state to get the initial pressure (if necessary)
    initial_pressure(cells[i]);

    // use the cell volume to convert primitive into conserved variables
    cells[i]._m = cells[i]._rho * cells[i]._V;
    cells[i]._p = cells[i]._m * cells[i]._u;
    cells[i]._E = cells[i]._P * cells[i]._V / (GAMMA - 1.) +
                  0.5 * cells[i]._u * cells[i]._p;
    const double cs =
        std::sqrt(GAMMA * cells[i]._P / cells[i]._rho) + std::abs(cells[i]._u);
    const double dt = courant_factor * cells[i]._V / cs;
    const uint_fast64_t integer_dt = (dt / maxtime) * integer_maxtime;
    min_integer_dt = std::min(min_integer_dt, integer_dt);

    cells[i]._last_rho = cells[i]._rho;
    cells[i]._last_u = cells[i]._u;
    cells[i]._last_P = cells[i]._P;
    cells[i]._last_nfac = cells[i]._nfac;
    cells[i]._last_entry = 0;
  }

  //  std::ofstream logfile("logfile.dat");
  LogFile logfile;
  initialize_log(logfile, "logfile.dat", MB_to_bytes(100));
  write_logfile(logfile, cells, ncell, 0., true);

  // set cell time steps
  // round min_integer_dt to closest smaller power of 2
  uint_fast64_t global_integer_dt = round_power2_down(min_integer_dt);
#pragma omp parallel for
  for (uint_fast32_t i = 1; i < ncell + 1; ++i) {
    cells[i]._integer_dt = global_integer_dt;
    cells[i]._dt = cells[i]._integer_dt * time_conversion_factor;
  }

  // initialize boundary condition and ionization variables
  boundary_conditions_initialize();
  ionization_initialize();

  // initialize the Riemann solver
  HLLCRiemannSolver solver(GAMMA);

  Timer step_time;
  double time_since_last = 0.;
  double time_since_start = 0.;
  unsigned int steps_since_last = 0;
  uint_fast64_t current_integer_time = 0;
  uint_fast64_t current_integer_dt = global_integer_dt;
  // main simulation loop: perform NSTEP steps
  const uint_fast64_t snaptime = DT * SNAPSTEP / maxtime * integer_maxtime;
  uint_fast64_t isnap = 0;
  while (current_integer_time < integer_maxtime) {

    step_time.start();

    add_spherical_source_term();

    // do first gravity kick
    do_gravity();

    // do ionization
    do_ionization();

    // update the primitive variables based on the values of the conserved
    // variables and the current cell volume
    min_integer_dt = integer_maxtime;
#pragma omp parallel for reduction(min : min_integer_dt)
    for (uint_fast32_t i = 1; i < ncell + 1; ++i) {
      cells[i]._rho = cells[i]._m / cells[i]._V;
      cells[i]._u = cells[i]._p / cells[i]._m;
      update_pressure(cells[i]);
      const double cs = std::sqrt(GAMMA * cells[i]._P / cells[i]._rho) +
                        std::abs(cells[i]._u);
      const double dt = courant_factor * cells[i]._V / cs;
      const uint_fast64_t integer_dt = (dt / maxtime) * integer_maxtime;
      min_integer_dt = std::min(min_integer_dt, integer_dt);
    }

    write_logfile(logfile, cells, ncell,
                  current_integer_time * time_conversion_factor);

    // round min_integer_dt to closest smaller power of 2
    global_integer_dt = round_power2_down(min_integer_dt);
    // make sure the remaining time can be filled *exactly* with the current
    // time step
    while ((integer_maxtime - current_integer_time) % global_integer_dt > 0) {
      global_integer_dt >>= 1;
    }
#pragma omp parallel for
    for (uint_fast32_t i = 1; i < ncell + 1; ++i) {
      cells[i]._integer_dt = global_integer_dt;
      cells[i]._dt = cells[i]._integer_dt * time_conversion_factor;
    }
    current_integer_dt = global_integer_dt;

    // check if we need to output a snapshot
    if (current_integer_time >= isnap * snaptime) {
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
      std::cout << "\t\t\tCentral mass: " << central_mass << " ("
                << (central_mass / MASS_POINT_MASS) << ")" << std::endl;
#endif
      time_since_last = 0.;
      steps_since_last = 0;
      write_snapshot(isnap, current_integer_time * time_conversion_factor,
                     cells, ncell);
      ++isnap;
    }

    // apply boundary conditions
    // SHOULD BE FAST ENOUGH NOT TO CARE ABOUT PARALLELIZATION
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
          std::min(1., 0.5 * std::min((rhomax - cells[i]._rho) / rhoextmax,
                                      (rhomin - cells[i]._rho) / rhoextmin));
      cells[i]._grad_rho = alpha_rho * gradrho;

      const double gradu = (cells[i + 1]._u - cells[i - 1]._u) * dx_inv;
      const double umax = std::max(cells[i - 1]._u, cells[i + 1]._u);
      const double umin = std::min(cells[i - 1]._u, cells[i + 1]._u);
      const double u_ext_plu = half_dx * gradu;
      const double u_ext_min = -half_dx * gradu;
      const double uextmax = std::max(u_ext_min, u_ext_plu);
      const double uextmin = std::min(u_ext_min, u_ext_plu);
      const double alpha_u =
          std::min(1., 0.5 * std::min((umax - cells[i]._u) / uextmax,
                                      (umin - cells[i]._u) / uextmin));
      cells[i]._grad_u = alpha_u * gradu;

      const double gradP = (cells[i + 1]._P - cells[i - 1]._P) * dx_inv;
      const double Pmax = std::max(cells[i - 1]._P, cells[i + 1]._P);
      const double Pmin = std::min(cells[i - 1]._P, cells[i + 1]._P);
      const double P_ext_plu = half_dx * gradP;
      const double P_ext_min = -half_dx * gradP;
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
      cells[i]._u -= half_dt * (u * cells[i]._grad_u + cells[i]._grad_P / rho);
      cells[i]._P -=
          half_dt * (GAMMA * P * cells[i]._grad_u + u * cells[i]._grad_P);

      if (cells[i]._rho < 0.) {
        cells[i]._rho = rho;
      }
      if (cells[i]._P < 0.) {
        cells[i]._P = P;
      }

      add_gravitational_prediction(cells[i], half_dt);
    }

// do the flux exchange
#pragma omp parallel for
    for (uint_fast32_t i = 1; i < ncell + 1; ++i) {
      const double dt = cells[i]._dt;
      // left flux
      {
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

        double mflux, pflux, Eflux;
        solver.solve_for_flux(rhoL_dash, uL_dash, PL_dash, rhoR_dash, uR_dash,
                              PR_dash, mflux, pflux, Eflux);

        cells[i]._m += dt * mflux;
        cells[i]._p += dt * pflux;
        cells[i]._E += dt * Eflux;
        if (i == 1) {
          flux_into_inner_mask(dt * mflux);
        }
      }
      // right flux
      {
        // get the variables in the left and right state
        const double rhoL = cells[i]._rho;
        const double uL = cells[i]._u;
        const double PL = cells[i]._P;
        const double rhoR = cells[i + 1]._rho;
        const double uR = cells[i + 1]._u;
        const double PR = cells[i + 1]._P;

        // do the second order spatial reconstruction
        const double dmin = 0.5 * (cells[i + 1]._midpoint - cells[i]._midpoint);
        const double dplu = -dmin;
        double rhoL_dash, rhoR_dash, uL_dash, uR_dash, PL_dash, PR_dash;
        rhoL_dash = rhoL + dmin * cells[i]._grad_rho;
        uL_dash = uL + dmin * cells[i]._grad_u;
        PL_dash = PL + dmin * cells[i]._grad_P;
        rhoR_dash = rhoR + dplu * cells[i + 1]._grad_rho;
        uR_dash = uR + dplu * cells[i + 1]._grad_u;
        PR_dash = PR + dplu * cells[i + 1]._grad_P;

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

        double mflux, pflux, Eflux;
        solver.solve_for_flux(rhoL_dash, uL_dash, PL_dash, rhoR_dash, uR_dash,
                              PR_dash, mflux, pflux, Eflux);

        cells[i]._m -= dt * mflux;
        cells[i]._p -= dt * pflux;
        cells[i]._E -= dt * Eflux;
      }
    }

    add_spherical_source_term();

    // do second gravity kick
    do_gravity();

    step_time.stop();
    time_since_last += step_time.value();
    ++steps_since_last;
    step_time.reset();
    current_integer_time += current_integer_dt;
  }

  // write the final logfile entry
  write_logfile(logfile, cells, ncell,
                current_integer_time * time_conversion_factor, true);

  close_log(logfile);

  // write the final snapshot
  write_snapshot(isnap, current_integer_time * time_conversion_factor, cells,
                 ncell);
  write_binary_snapshot(cells, ncell);

  delete[] cells;

  total_time.stop();
  std::cout << "Total program time: " << total_time.value() << " s."
            << std::endl;

  // all went well: return with exit code 0
  return 0;
}
