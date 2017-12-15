#! /bin/bash

# run a Sod shock test in the folder set below, using the number of threads set
# below

folder=build_sod
nthread=8

cmake_command=$(python write_configuration_sod.py)

mkdir $folder
cd $folder
echo $cmake_command
eval $cmake_command
make -j $nthread
OMP_NUM_THREADS=$nthread OMP_PROC_BIND=True ./HydroCodeSpherical1D
cd ..
