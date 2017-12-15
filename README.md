# HydroCodeSpherical1D

1D spherical hydrodynamics code with optional external point mass and two
temperature approach ionisation.

To run the code, download a copy of this repository. Then, configure the code
using `cmake` in a folder of choice:
```
cmake -DCMAKE_BUILD_TYPE=Release PATH_TO_SOURCE
```
where `PATH_TO_SOURCE` is the path to the folder containing `CMakeLists.txt`.
We strongly recommend configuring the code in a different folder from
`PATH_TO_SOURCE` to avoid problems when configuring different versions. An
example configuration can also be obtained by running either `run_sod_shock.sh`,
`run_bondi_setup.sh` or `run_bondi_run.sh`, which configures, compiles and runs
the code in a subfolder of the main project folder. To help setting up a
configuration, we provided the script `get_cmake_command.py`, which can also be
used as a Python module, see e.g. `write_configuration_sod.py`.

Once configured, run `make` from the folder where `cmake` was run:
```
make
```
This will compile the program. Note that the program consists of a single file,
so trying to speed up the compilation (using `make -j 4` to e.g. build with 4
threads) will not work.

Once successfully compiled, the program can be run using
```
./HydroCodeSpherical1D
```
from the same folder where `make` and `cmake` were run. The program will provide
ample output during the run, and will, depending on the configuration, create a
(large) number of output files named `snapshot_XXXX.txt`. These contain the
midpoint position, density, velocity, pressure and neutral fraction for each
cell in the grid (in SI units). The program will also create a file named
`lastsnap.dat`, which can be used as initial condition for a new run (if the
code is configured with `ic_mode=IC_FILE`). If self-consistent ionisation is
used (configuration option `ionisation_mode=IONISATION_MODE_SELF_CONSISTENT`),
an output file `ionisation_radius.dat` is also created, which contains a binary
dump of the ionisation radius as a function of time (using a smart conditional
output algorithm). Finally, an empty file named `logfile.dat` is also created.
This file can be safely ignored.
