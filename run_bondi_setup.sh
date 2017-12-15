#! /bin/bash

# run a Bondi setup simulation in the folder given below, with the number of
# threads set below

folder=build_bondi_setup
nthread=8

cmake_command=$(python write_configuration_bondi_setup.py)

mkdir $folder
cd $folder
echo $cmake_command
eval $cmake_command
make -j $nthread
OMP_NUM_THREADS=$nthread OMP_PROC_BIND=True ./HydroCodeSpherical1D 2>&1 \
  | tee bondi_setup.log
cd ..
