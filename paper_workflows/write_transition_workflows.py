import sys
sys.path.append("..")
import get_cmake_command

git_command = \
  "git clone https://github.com/bwvdnbro/HydroCodeSpherical1D.git source"

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

make_command = "make"
run_command = \
  "OMP_NUM_THREADS=8 OMP_PROC_BIND=True ./HydroCodeSpherical1D > run.log 2>&1"
icname = "convergence_stable_w{width:.0f}_{resolution}.dat"
snapname = "convergence_stable_w{width:.0f}_2700.txt"
instable_radius_name = \
  "convergence_instability_w{width:.0f}_{resolution}_radius.dat"
seed_radius_name = \
  "convergence_seed_w{width:.0f}_{resolution}_{sign}{amplitude}_radius.dat"
seed_command = \
  "python source/add_density_perturbation.py ic_noseed.dat ic.dat {amplitude}"

for width in [1., 2., 3., 4., 5.]:
  for resolution in [300, 900, 2700, 5400]:
    bondi_setup_this = dict(bondi_setup)
    bondi_setup_this["ncell"] = resolution
    bondi_setup_this["ionisation_transition_width_in_au"] = width
    cmake_setup = get_cmake_command.get_cmake_command(bondi_setup_this,
                                                      "source/")
    icname_this = icname.format(width = width, resolution = resolution)
    output = "{0}->lastsnap.dat".format(icname_this)
    if resolution == 2700:
      snapname_this = snapname.format(width = width)
      output += " {0}->snapshot_2000.txt".format(snapname_this)

    print "{output}:".format(output = output)
    print "\t{git}; {cmake}; {make}; {run}\n".format(
      git = git_command, cmake = cmake_setup, make = make_command,
      run = run_command)

    bondi_run_this = dict(bondi_run)
    bondi_run_this["ncell"] = resolution
    bondi_run_this["ionisation_transition_width_in_au"] = width
    cmake_run = get_cmake_command.get_cmake_command(bondi_run_this, "source/")
    instable_radius_name_this = instable_radius_name.format(
      width = width, resolution = resolution)
    output = "{0}->ionisation_radius.dat".format(instable_radius_name_this)
    input = "{0}->ic.dat".format(icname_this)
    print "{output}: {input}".format(output = output, input = input)
    print "\t{git}; {cmake}; {make}; {run}\n".format(
      git = git_command, cmake = cmake_run, make = make_command,
      run = run_command)

    for amplitude in [0.1, 2., -0.1]:
      seed_radius_name_this = seed_radius_name.format(
        width = width, resolution = resolution,
        sign = 'p' if amplitude > 0. else 'm', amplitude = abs(amplitude))
      output = "{0}->ionisation_radius.dat".format(seed_radius_name_this)
      if resolution == 2700 and width == 5.0:
        for snap in [0, 325, 425, 650, 900]:
          output += \
            " convergence_seed_{0:03d}_{1}{2}.txt->snapshot_{0:04d}.txt".format(
              snap, 'p' if amplitude > 0. else 'm', abs(amplitude))
      input = "{0}->ic_noseed.dat".format(icname_this)
      seed_command_this = seed_command.format(amplitude = amplitude)
      print "{output}: {input}".format(output = output, input = input)
      print "\t{git}; {cmake}; {make}; {seed}; {run}\n".format(
        git = git_command, cmake = cmake_run, make = make_command,
        seed = seed_command_this, run = run_command)
