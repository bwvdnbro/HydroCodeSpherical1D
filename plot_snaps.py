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
# @file plot_snaps.py
#
# @brief Reads all snapshot files in the current folder and makes a plot of the
# density, velocity, pressure and neutral fraction profile as a function of
# radius.
#
# @autor Bert Vandenbroucke (bv7@st-andrews.ac.uk)
##

import numpy as np
import pylab as pl
import glob
import sys
import multiprocessing as mp

# unit conversions: we plot distances in AU and give time in years
au_in_si = 1.495978707e11 # m
yr_in_si = 365.25 * 24. * 3600. # s

# number of parallel threads to use: given as a command line parameter
# (default: 1)
nthread = 1
if len(sys.argv) > 1:
  nthread = int(sys.argv[1])

##
# @brief Make the plots for a single snapshot file.
#
# @param f Snapshot file to plot.
# @return Name of the file.
##
def plot(f):
  # read the time stamp at the start of the snapshot file
  ifile = open(f, 'r')
  timeline = ifile.readline()
  time = float(timeline.split()[2])
  ifile.close()
  # now read in the actual data values
  data = np.loadtxt(f)
  # unit conversion
  data[:,0] /= au_in_si

  # make the plots
  fig, ax = pl.subplots(2, 2, sharex = "col")
  ax[0][0].plot(data[:,0], data[:,1]*0.001, "k-")
  ax[0][1].plot(data[:,0], data[:,2]*0.001, "k-")
  ax[1][0].plot(data[:,0], data[:,3], "k-")
  ax[1][1].plot(data[:,0], data[:,4], "k-")
  ax[0][0].set_ylabel("density (g cm$^{-3}$)")
  ax[0][1].set_ylabel("velocity (km s$^{-1}$)")
  ax[1][0].set_ylabel("pressure (kg m$^{-1}$ s$^{-2}$)")
  ax[1][1].set_ylabel("neutral fraction")
  ax[1][0].set_xlabel("radius (AU)")
  ax[1][1].set_xlabel("radius (AU)")
  ax[0][1].set_title("t = {t:.2e} yr".format(t = time / yr_in_si))
  pl.tight_layout()
  pl.savefig("{name}.png".format(name = f[:-4]))
  pl.close()

  return f

# set up a parallel pool to do the plotting
pool = mp.Pool(nthread)
results = []
# scan the working directory for snapshot files
for f in sorted(glob.glob("snapshot_*.txt")):
  # add the file to the list of tasks
  results.append(pool.apply_async(plot, (f,)))

# wait for the threads to finish
for result in results:
  f = result.get()
  print "Done plotting", f
