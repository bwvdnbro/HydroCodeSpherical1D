#! /bin/bash

# run a Sedov blastwave test in the folder set below, using the number of
# threads set below

folder=build_sedov_cooling
nthread=1

cmake_command=$(python3 write_configuration_sedov_cooling.py)

mkdir $folder
cd $folder
echo $cmake_command
eval $cmake_command
make -j $nthread
OMP_NUM_THREADS=$nthread OMP_PROC_BIND=True ./HydroCodeSpherical1D
cd ..
