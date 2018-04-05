#! /usr/bin/python

################################################################################
# This file is part of HydroCodeSpherical1D
# Copyright (C) 2017, 2018 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
# @file plot_starbench.py
#
# @brief Plot the binary ionisation radius as a function of time file.
#
# @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
##

import numpy as np
import matplotlib
matplotlib.use("Agg")
import pylab as pl

# units: we plot time in Myr and distance in pc
pc_in_si = 3.086e16
Myr_in_si = (1.e6 * 365.25 * 24. * 3600.)

# initial ionisation radius
au_in_si = 1.496e11
Rst_in_si = 6.5e4 * au_in_si

# ionised sound speed
ci_in_si = 1.285e4

# time range for reference solutions
t = np.linspace(0., 0.141, 1000) * Myr_in_si

# Spitzer solution
Rsp_in_si = Rst_in_si * (1. + 1.75 * ci_in_si * t / Rst_in_si)**(4. / 7.)

# Hosokawa-Inutsuka solution
Rhi_in_si = Rst_in_si * \
            (1. + 1.75 * np.sqrt(4. / 3.) * ci_in_si * t / Rst_in_si)**(4. / 7.)

file = "ionisation_radius.dat"

# memory-map the binary file to a read-only numpy array
fp = np.memmap(file, dtype = 'd', mode = 'r')
# the file has 3 columns: the time, ionisation radius and ionising luminosity
# ratio (in SI units)
data = fp.reshape((-1, 3))

# create the plot
fig, ax = pl.subplots(1, 1, sharex = True, sharey = True)
ax.plot(t / Myr_in_si, Rsp_in_si / pc_in_si, "r-", label = "Spitzer")
ax.plot(t / Myr_in_si, Rhi_in_si / pc_in_si, "b-", label = "Hosokawa-Inutsuka")
ax.plot(data[:,0] / Myr_in_si, data[:,1] / pc_in_si, "k.", label = "simulation")
ax.set_ylabel("ionisation radius (pc)")
ax.set_xlabel("time (Myr)")
ax.legend(loc = "best")
pl.tight_layout()
# show the plot window
pl.savefig("starbench.png")
