#! /usr/bin/python

################################################################################
# This file is part of HydroCodeSpherical1D
# Copyright (C) 2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
#
# HydroCodeSpherical1D is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# HydroCodeSpherical1D is distributed in the hope that it will be useful,
# but WITOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with HydroCodeSpherical1D. If not, see <http://www.gnu.org/licenses/>.
################################################################################

##
# @file add_density_perturbation.py
#
# @brief Script that adds a density perturbation to the given binary output file
# and stores the result in a new binary output file with the given name.
#
# @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
##

import numpy as np
import sys

##
# @brief Cubic spline kernel like bump filter.
#
# Based on the cubic spline expression in Springel (2005).
#
# @param x Position array.
# @param center Central position of the bump.
# @param h Width of the bump.
# @param factor Amplitude of the bump.
# @return Filter array that is 1 for every x value outside the bump region, and
# 1 + bump inside the bump region.
##
def cubic_spline(x, center, h, factor):
  return 1. + factor * \
         np.where(abs(x - center) < 0.5 * h,
                  1. - 6. * ((x - center) / h)**2 + \
                    6. * (abs(x - center) / h)**3,
                  np.where(abs(x - center) < h,
                           2. * (1. - abs(x - center) / h)**3,
                           np.zeros(len(x))))                  

# Make sure we have 2 command line arguments: the input and output file name
if len(sys.argv) < 3:
  print "Usage: python add_density_perturbation input_file output_file [neg]"
  exit()

fac = 0.1

name = sys.argv[1]
outname = sys.argv[2]
if sys.argv[3] == "neg":
  fac = -fac

print fac

# memory-map the input file to a read-only numpy array
fp = np.memmap(name, dtype = 'd', mode = 'r')
# we know the array has 4 columns
data = fp.reshape((-1, 4))

# set up the spatial grid with the right resolution
r = np.linspace(10., 100., len(data) + 1)
r = 0.5 * (r[1:] + r[:-1])

# make a copy of the original data, which we will modify below
copy = data.copy()

# apply the density bump filter
copy[:,0] *= np.where(r > 60., np.where(r < 70., cubic_spline(r, 65., 5., fac),
                                                 np.ones(len(r))),
                               np.ones(len(r)))

# store the result in a new memory-mapped file
ofp = np.memmap(outname, dtype = 'd', mode = 'w+', shape = copy.shape)
ofp[:] = copy[:]
