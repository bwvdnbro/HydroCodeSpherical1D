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
# @file plot_ionisation_radius.py
#
# @brief Plot the binary ionisation radius as a function of time file.
#
# @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
##

import numpy as np
import pylab as pl

# units: we plot time in years and distance in AU
au_in_si = 1.496e11
yr_in_si = 365.25 * 24.0 * 3600.0

file = "ionisation_radius.dat"

# memory-map the binary file to a read-only numpy array
fp = np.memmap(file, dtype="d", mode="r")
# the file has 3 columns: the time, ionisation radius and ionising luminosity
# ratio (in SI units)
data = fp.reshape((-1, 3))

# create the plot
fig, ax = pl.subplots(1, 1, sharex=True, sharey=True)
ax.plot(data[:, 0] / yr_in_si, data[:, 1] / au_in_si)
ax.set_ylabel("ionisation radius (AU)")
ax.set_xlabel("time (yr)")
pl.tight_layout()
# show the plot window
pl.show()
