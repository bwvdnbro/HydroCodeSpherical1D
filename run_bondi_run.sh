#! /bin/bash

# run a Bondi stability test in the folder set below, with the number of threads
# set below

folder=build_bondi_run
nthread=8

# check if the initial condition (output of a Bondi setup run) exists
# if not, generate it by running a Bondi setup run first
if [ ! -e build_bondi_setup/lastsnap.dat ]
then
  ./run_bondi_setup.sh
fi

cmake_command=$(python write_configuration_bondi_run.py)

mkdir $folder
# add a well defined density perturbation to seed the instability
python add_density_perturbation.py build_bondi_setup/lastsnap.dat \
  $folder/bondi_ic.dat
cd $folder
echo $cmake_command
eval $cmake_command
make -j $nthread
OMP_NUM_THREADS=$nthread OMP_PROC_BIND=True ./HydroCodeSpherical1D 2>&1 \
  | tee bondi_run.log
cd ..
