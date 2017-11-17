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
 * @file LogFile.hpp
 *
 * @brief Memory-mapped log file output, as in SWIFT.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef LOGFILE_HPP
#define LOGFILE_HPP

#include <cstring>
#include <fcntl.h>
#include <iostream>
#include <string>
#include <sys/mman.h>
#include <unistd.h>

/**
 * @brief Memory-mapped log file.
 */
class LogFile {
private:
  /*! @brief The memory-mapped part of the log file. */
  char *_memory_buffer;

  /*! @brief The number of bytes of the memory-mapped part of the file that have
   *  already been used for output, in bytes. */
  size_t _memory_buffer_count;

  /*! @brief The offset of the memory-mapped data within the log file, in
   *  bytes. */
  size_t _file_offset;

  /*! @brief The log file. */
  int _file;

  /*! @brief Mask used to ensure sizes and offsets that are multiples of the
   *  system page size. */
  const size_t _page_mask;

  /*! @brief The size of the memory-mapped data, in bytes. */
  const size_t _memory_buffer_size;

  /**
   * @brief Get the page size mask that can be used to round sizes and offsets
   * to page size multiples.
   *
   * To round a variable `s` up, use
   * ```
   *  s = (s + ~page_mask) & page_mask;
   * ```
   *
   * To round down, use
   * ```
   *  s &= page_mask;
   * ```
   *
   * Examples (we assume a 16-bit mask for simplicity, page size is \f$PS =
   * 4096\f$, page mask is \f$PM = 1111~0000~0000~0000\f$):
   *  - \f$s = 5 = 101\f$:
   *      \f[
   *        (s + \sim{}PM) \& PM
   *          = (101 + 0000~1111~1111~1111) & 1111~0000~0000~0000
   *          = 0001~0000~0000~0100 & 1111~0000~0000~0000
   *          = 0001~0000~0000~0000
   *      \f]
   *  - \f$s = PM + 2\f$:
   *      \f[
   *        (s + \sim{}PM) \& PM
   *          = (10 + 0001~0000~0000~0000 + 0000~1111~1111~1111) &
   *            1111~0000~0000~0000
   *          = 0010~0000~0000~0100 & 1111~0000~0000~0000
   *          = 0010~0000~0000~0000
   *      \f]
   *
   * @return Page mask that can be used to round sizes and offsets to page size
   * multiples.
   */
  static inline size_t get_page_mask() { return ~(sysconf(_SC_PAGE_SIZE) - 1); }

  /**
   * @brief Round the given size up to a multiple of the page size.
   *
   * @param size Size to round up.
   */
  inline size_t round_page_up(const size_t size) const {
    return (size + ~_page_mask) & _page_mask;
  }

  /**
   * @brief Round the given size down to a multiple of the page size.
   *
   * @param size Size to round down.
   */
  inline size_t round_page_down(const size_t size) const {
    return size & _page_mask;
  }

  /**
   * @brief Increase the log file size by shifting the part of the file that was
   * written to disk.
   *
   * We unmap the part of the file that was already written and grow the file to
   * that size plus the default buffer size. We then memory-map the new part of
   * the file.
   */
  inline void increase_file_size() {
    // truncate the number of bytes already written to the file to a multiple of
    // the page size
    // we will shift the memory-mapped region to this point
    const size_t buffer_offset = round_page_down(_memory_buffer_count);
    // unmap the part of the memory-mapped region before the new offset
    // if nothing was written to the file, we provide 1, which is rounded up to
    // a single page size
    if (munmap(_memory_buffer, buffer_offset > 0 ? buffer_offset : 1) != 0) {
      std::cerr << "Error unmapping part of log file!" << std::endl;
      abort();
    }
    // increase the file offset, which indicates the offset of the memory-mapped
    // region within the file on disk
    _file_offset += buffer_offset;
    // subtract the unmapped region from the counter that counts the bytes
    // already written to the mapped region
    _memory_buffer_count -= buffer_offset;
    // make sure we have enough disk space for the next bit of log file
    if (posix_fallocate(_file, _file_offset, _memory_buffer_size) != 0) {
      std::cerr << "Error reserving extra log file size on disk!" << std::endl;
      abort();
    }
    // map the new bit of file
    // parameters are the same as above, but now we start from a later offset in
    // the file
    _memory_buffer =
        reinterpret_cast<char *>(mmap(nullptr, _memory_buffer_size, PROT_WRITE,
                                      MAP_SHARED, _file, _file_offset));
    if (_memory_buffer == MAP_FAILED) {
      std::cerr << "Error memory mapping new part of log file!" << std::endl;
      abort();
    }
  }

public:
  /**
   * @brief Constructor.
   *
   * Creates a log file with the given file name and size, and memory-maps it
   * for fast write access.
   *
   * @param filename Name of the file.
   * @param size Initial size of the file, in MB. The same size will be used for
   * the internal file buffer that is stored in memory.
   */
  inline LogFile(const std::string filename, const size_t size)
      : _memory_buffer(nullptr), _memory_buffer_count(0), _file_offset(0),
        _file(0), _page_mask(get_page_mask()),
        _memory_buffer_size(round_page_up(size << 20)) {
    // create the file
    // O_CREAT: create the file if it does not exist
    // O_RDWR: read/write access
    // S_IRUSR/IWUSR/IRGRP/IWGRP: set read/write access for both the user and
    //  group when creating the file (code 0660 in the SWIFT code)
    _file = open(filename.c_str(), O_CREAT | O_RDWR,
                 S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP);
    if (_file < 0) {
      std::cerr << "Error opening log file!" << std::endl;
      abort();
    }
    // make sure we have enough space on disk to store the memory mapped file
    // this preallocates the file size without actually writing anything
    // we tell the call we want all file offsets between 0 and size to be
    // available after this call
    if (posix_fallocate(_file, 0, _memory_buffer_size) != 0) {
      std::cerr << "Error reserving log file disk space!" << std::endl;
      abort();
    }
    // memory-map the file: make a synced copy of some part of the file in
    // active memory that we can directly access as if it was normal memory
    // parameters:
    //  address: let the system decide where to allocate the memory
    //  size: size of the buffer we want to memory-map (should be a multiple of
    //   page_size!)
    //  protection: we want to write to the buffer
    //  mode: we want to share the memory with all other processes: every write
    //   to this region is synced in the file and directly reflected in all
    //   other processes that memory map the same region
    //  file: this is the file we want to memory-map
    //  offset: we want to start memory mapping from offset 0 in the file
    _memory_buffer = reinterpret_cast<char *>(
        mmap(nullptr, _memory_buffer_size, PROT_WRITE, MAP_SHARED, _file, 0));
    if (_memory_buffer == MAP_FAILED) {
      std::cerr << "Error memory mapping log file!" << std::endl;
      abort();
    }
  }

  /**
   * @brief Close the log file.
   *
   * Unmaps the memory-mapped part of the file (forcing it to be written to the
   * file) and closes the log file.
   */
  inline void close_file() {
    // unmap the actively memory-mapped region of the file
    if (munmap(_memory_buffer,
               _memory_buffer_count > 0 ? _memory_buffer_count : 1) != 0) {
      std::cerr << "Error unmapping log file memory!" << std::endl;
      abort();
    }
    _memory_buffer = nullptr;
    // shrink the file to its actual size
    if (ftruncate(_file, _file_offset + _memory_buffer_count) != 0) {
      std::cerr << "Error shrinking log file!" << std::endl;
      abort();
    }
    // close the file
    if (close(_file) != 0) {
      std::cerr << "Error closing log file!" << std::endl;
      abort();
    }
    _file = 0;
  }

  /**
   * @brief Destructor.
   */
  inline ~LogFile() {
    if (_file != 0) {
      close_file();
    }
  }

  /**
   * @brief Make sure the file has enough space to contain the given number of
   * bytes and grow the file and shift the memory buffer if necessary.
   *
   * @param size Minimum requested available size, in bytes.
   */
  inline void ensure_space(const size_t size) {
    if (_memory_buffer_size < _memory_buffer_count + size) {
      increase_file_size();
    }
  }

  /**
   * @brief Get the current writing position in the file.
   *
   * @return Current writing position in the file, in bytes.
   */
  inline size_t get_current_position() const {
    return _file_offset + _memory_buffer_count;
  }

  /**
   * @brief Make sure all contents in the log file are actually written to disk.
   */
  inline void flush() {
    // make sure the memory-mapped region is synced on disk
    // parameters:
    //  address: pointer to the beginning of the region we want synced
    //  size: number of bytes after that pointer we want to be synced
    //  MS_SYNC: force actual writing to disk and wait until it is finished
    msync(_memory_buffer, _memory_buffer_count, MS_SYNC);
  }

  /**
   * @brief Write the given value to the log file.
   *
   * @param value Value to write.
   */
  template <typename _type_> inline void write(_type_ value) {
    const size_t vsize = sizeof(_type_);
    ensure_space(vsize);
    const size_t offset = _memory_buffer_count;
    _memory_buffer_count += vsize;
    memcpy(_memory_buffer + offset, &value, vsize);
  }
};

#endif // LOGFILE_HPP
