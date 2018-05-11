import sys
sys.path.append("..")
import get_cmake_command

# command required to get a fresh copy of the repository from GitHub
git_command = \
  "git clone https://github.com/bwvdnbro/HydroCodeSpherical1D.git source"

# code configuration for a Bondi setup run: homogeneous initial density with
# constant inflow and a fixed ionisation radius
bondi_setup = {
"rmin_in_au": 10.,
"rmax_in_au": 100.,
"ncell": 2700,
"gamma": 1.001,
"maxtime_in_yr": 40.,
"number_of_snaps": 2000,
"ic": "IC_BONDI",
"eos": "EOS_BONDI",
"boundaries": "BOUNDARIES_BONDI",
"isothermal_temperature_in_k": 500.,
"potential": "POTENTIAL_POINT_MASS",
"g_internal": 1.,
"mass_point_mass_in_msol": 18.,
"bondi_density_in_si": 1.e-16,
"bondi_pressure_contrast": 32.,
"initial_ionisation_radius_in_au": 30.,
"unit_mass_in_si": 2.479e31,
"unit_length_in_si": 1.2e13,
"ionisation_mode": "IONISATION_MODE_CONSTANT",
"ionisation_transition": "IONISATION_TRANSITION_SMOOTH",
"ionisation_transition_width_in_au": 5.,
"courant_factor": 0.05,
"riemannsolver_type": "RIEMANNSOLVER_TYPE_HLLC"
}

# code configuration for a Bondi stability run: density field read from a file
# with self-consistent ionisation
bondi_run = {
"rmin_in_au": 10.,
"rmax_in_au": 100.,
"ncell": 2700,
"gamma": 1.001,
"maxtime_in_yr": 80.,
"number_of_snaps": 4000,
"ic": "IC_FILE",
"eos": "EOS_BONDI",
"boundaries": "BOUNDARIES_BONDI",
"isothermal_temperature_in_k": 500.,
"potential": "POTENTIAL_POINT_MASS",
"g_internal": 1.,
"mass_point_mass_in_msol": 18.,
"bondi_density_in_si": 1.e-16,
"bondi_pressure_contrast": 32.,
"ic_file_name": "ic.dat",
"initial_ionisation_radius_in_au": 30.,
"unit_mass_in_si": 2.479e31,
"unit_length_in_si": 1.2e13,
"ionisation_mode": "IONISATION_MODE_SELF_CONSISTENT",
"ionisation_transition": "IONISATION_TRANSITION_SMOOTH",
"ionisation_transition_width_in_au": 5.,
"courant_factor": 0.05,
"riemannsolver_type": "RIEMANNSOLVER_TYPE_HLLC"
}

# compilation command
make_command = "make"
# run command for a simulation
run_command = \
  "OMP_NUM_THREADS=8 OMP_PROC_BIND=True ./HydroCodeSpherical1D > run.log 2>&1"
# initial condition file name for a run with a smooth ionisation transition
icname = "convergence_stable_w{width:.0f}_{resolution}.dat"
# final snapshot name for a representative IC run with a smooth transition
snapname = "convergence_stable_w{width:.0f}_2700.txt"
# ionisation radius log file name for a run with numerically seeded instability
instable_radius_name = \
  "convergence_instability_w{width:.0f}_{resolution}_radius.dat"
# ionisation radius log file name for a run with seeded instability
seed_radius_name = \
  "convergence_seed_w{width:.0f}_{resolution}_{sign}{amplitude}_radius.dat"
# command necessary to seed an instability
seed_command = \
  "python source/add_density_perturbation.py ic_noseed.dat ic.dat {amplitude}"

# workflow for all simulations that include a smooth transition and seeded
# instabilities
seed_lines = ""
# loop over all transition widths
for width in [1., 2., 3., 4., 5.]:
  # loop over all resolutions
  for resolution in [300, 900, 2700, 5400]:
  
    ## IC generation
  
    # make a copy of the Bondi setup configuration and set the right resolution
    # and transition width
    bondi_setup_this = dict(bondi_setup)
    bondi_setup_this["ncell"] = resolution
    bondi_setup_this["ionisation_transition_width_in_au"] = width
    # get the corresponding configuration command
    cmake_setup = get_cmake_command.get_cmake_command(bondi_setup_this,
                                                      "source/")
    # get the appropriate name for the final result (this will serve as IC for
    # the simulations below)
    icname_this = icname.format(width = width, resolution = resolution)
    output = "{0}->lastsnap.dat".format(icname_this)
    # for the 2700 cell runs: save the last snapshot for the profile plot
    if resolution == 2700:
      snapname_this = snapname.format(width = width)
      output += " {0}->snapshot_2000.txt".format(snapname_this)

    # add the IC task to the workflow
    seed_lines += "{output}:\n".format(output = output)
    seed_lines +=  "\t{git}; {cmake}; {make}; {run}\n\n".format(
      git = git_command, cmake = cmake_setup, make = make_command,
      run = run_command)

    ## stability test: test how long the transition suppresses numerical
    ## instabilities

    # make a copy of the Bondi run configuration and set the right resolution
    # and transition width
    bondi_run_this = dict(bondi_run)
    bondi_run_this["ncell"] = resolution
    bondi_run_this["ionisation_transition_width_in_au"] = width
    # get the corresponding configuration command
    cmake_run = get_cmake_command.get_cmake_command(bondi_run_this, "source/")
    # get the appropriate name for the ionisation radius log file
    instable_radius_name_this = instable_radius_name.format(
      width = width, resolution = resolution)
    output = "{0}->ionisation_radius.dat".format(instable_radius_name_this)
    input = "{0}->ic.dat".format(icname_this)

    # add the stability run task to the workflow
    seed_lines += "{output}: {input}\n".format(output = output, input = input)
    seed_lines += "\t{git}; {cmake}; {make}; {run}\n\n".format(
      git = git_command, cmake = cmake_run, make = make_command,
      run = run_command)

    ## seeded instability runs: seed an instability and analyze the behaviour

    # loop over all bump sizes
    for amplitude in [0.1, 2., -0.1, -0.01]:
      # only do a full width/resolution scan for the -0.1 runs
      if amplitude != -0.1 and (resolution != 2700 or width != 5.0):
        continue

      # get the appropriate name for the ionisation radius log file
      seed_radius_name_this = seed_radius_name.format(
        width = width, resolution = resolution,
        sign = 'p' if amplitude > 0. else 'm', amplitude = abs(amplitude))
      output = "{0}->ionisation_radius.dat".format(seed_radius_name_this)
      # add extra snapshots for the 2700 cell runs with a 5 AU transition
      if resolution == 2700 and width == 5.0:
        for snap in [0, 325, 425, 650, 900]:
          output += \
            " convergence_seed_{0:03d}_{1}{2}.txt->snapshot_{0:04d}.txt".format(
              snap, 'p' if amplitude > 0. else 'm', abs(amplitude))
      input = "{0}->ic_noseed.dat".format(icname_this)
      
      # add the simulation run task to the workflow
      seed_command_this = seed_command.format(amplitude = amplitude)
      seed_lines += "{output}: {input}\n".format(output = output, input = input)
      seed_lines += "\t{git}; {cmake}; {make}; {seed}; {run}\n\n".format(
        git = git_command, cmake = cmake_run, make = make_command,
        seed = seed_command_this, run = run_command)

print seed_lines
