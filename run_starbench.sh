#! /bin/bash

# run a starbench simulation in the folder given below, with the number of
# threads set below

folder=build_starbench
nthread=8

cmake_command=$(python write_configuration_starbench.py)

mkdir $folder
cd $folder
echo $cmake_command
eval $cmake_command
make -j $nthread
OMP_NUM_THREADS=$nthread OMP_PROC_BIND=True ./HydroCodeSpherical1D 2>&1 \
  | tee starbench.log
cd ..
