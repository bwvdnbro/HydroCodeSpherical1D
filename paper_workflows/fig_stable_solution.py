import numpy as np
import matplotlib
matplotlib.use("Agg")
import pylab as pl
import scipy.special.lambertw as lambertw

pl.rcParams["text.usetex"] = True
pl.rcParams["figure.figsize"] = (6, 8)
pl.rcParams["font.size"] = 14
pl.rcParams["axes.labelsize"] = 18

gamma = 5. / 3.

au_in_si = 1.495978707e11 # m
yr_in_si = 365.25 * 24. * 3600. # s

## Bondi
# input unit parameters
unit_length_in_si = 1.2e13
unit_mass_in_si = 2.479e31
G_in_si = 6.67408e-11
k_in_si = 1.38064852e-23
mH_in_si = 1.674e-27
solar_mass_in_si = 1.9891e30

# derived units
unit_time_in_si = np.sqrt(unit_length_in_si**3 / (unit_mass_in_si * G_in_si))
unit_density_in_si = unit_mass_in_si / (unit_length_in_si**3)
unit_velocity_in_si = unit_length_in_si / unit_time_in_si
unit_pressure_in_si = unit_mass_in_si / (unit_length_in_si * unit_time_in_si**2)

# input parameters
# physical
mass_point_mass = 18. * solar_mass_in_si
T_n = 500
pressure_contrast = 32.
bondi_rho_n = 1.e-16
r_ion = 30. * au_in_si
# practical
r_min = 10. * au_in_si
r_max = 100. * au_in_si

# derived parameters
cs2_n = T_n * k_in_si / mH_in_si
cs2_i = pressure_contrast * cs2_n
bondi_r_n = 0.5 * G_in_si * mass_point_mass / cs2_n
bondi_r_i = 0.5 * G_in_si * mass_point_mass / cs2_i
cs_n = np.sqrt(cs2_n)
cs_i = np.sqrt(cs2_i)

def neutral_bondi(r):
  global cs_n, cs2_n, bondi_r_n, bondi_rho_n
  u = bondi_r_n / r
  omega = -u**4 * np.exp(3. - 4. * u)
  v = np.where(r < bondi_r_n,
                 -cs_n * np.sqrt(-lambertw(omega, -1).real),
                 -cs_n * np.sqrt(-lambertw(omega, 0).real))
  rho = -bondi_rho_n * bondi_r_n**2 * cs_n / r**2 / v
  P = cs2_n * rho
  return rho, v, P, r * 0. + 1.

rho1, v1, _, _ = neutral_bondi(r_ion)

Gamma = 0.5 * (v1**2 + cs2_n - \
               np.sqrt((v1**2 + cs2_n)**2 - 4. * v1**2 * cs2_i)) / cs2_i

rho2 = Gamma * rho1
v2 = v1 / Gamma

def ionised_bondi(r):
  global r_ion, v2, cs_i, cs2_i, bondi_r_i, rho2
  omega = -(v2 / cs_i)**2 * (r_ion / r)**4 * \
          np.exp(4. * (bondi_r_i / r_ion - bondi_r_i / r) - v2**2 / cs_i**2)
  v = -cs_i * np.sqrt(-lambertw(omega, -1).real)
  rho = r_ion**2 * v2 * rho2 / r**2 / v
  P = cs2_i * rho
  return rho, v, P, rho * 0.

ra = np.linspace(r_min, r_max, 1000)
rhoa, va, Pa, na = np.where(ra < r_ion, ionised_bondi(ra), neutral_bondi(ra))
ra /= au_in_si
rhoa *= 0.001
va *= 0.001

def plot(f, ax):
  ifile = open(f, 'r')
  timeline = ifile.readline()
  time = float(timeline.split()[2]) / yr_in_si
  ifile.close()
  data = np.loadtxt(f)
  data[:,0] /= au_in_si
  data[:,1] *= 0.001
  data[:,2] *= 0.001
  ax[0].semilogy(data[:,0], data[:,1], "-",
                 label = "$t = {0:.0f}~{{\\rm{{}}yr}}$".format(time))
  ax[1].plot(data[:,0], data[:,2], "-")

fig, ax = pl.subplots(2, 1, sharex = True)

plot("stable_solution_t00.txt", ax)
plot("stable_solution_t05.txt", ax)
plot("stable_solution_t10.txt", ax)
plot("stable_solution_t40.txt", ax)

ax[0].legend(loc = "best")

ax[0].plot(ra, rhoa, "k--", linewidth = 0.8)
ax[1].plot(ra, va, "k--", linewidth = 0.8)
ax[0].set_ylabel("$\\rho{}$ (g cm$^{-3}$)")
ax[1].set_ylabel("$v$ (km s$^{-1}$)")

ax[1].set_xlabel("$r$ (AU)")

pl.tight_layout()
pl.savefig("figure_stable_solution.png")
pl.close()

# relative error

data = np.loadtxt("stable_solution_t40.txt")

rhoa, va, Pa, na = np.where(data[:,0] < r_ion, ionised_bondi(data[:,0]),
                                               neutral_bondi(data[:,0]))

# unit conversion
data[:,0] /= au_in_si

pl.rcParams["figure.figsize"] = (6, 8)

fig, ax = pl.subplots(2, 1, sharex = True)

ax[0].semilogy(data[:,0], abs(rhoa - data[:,1]) / abs(rhoa + data[:,1]))
ax[1].semilogy(data[:,0], abs(va - data[:,2]) / abs(va + data[:,2]))
ax[0].set_ylabel("$|\\rho{}_s - \\rho{}_a| / |\\rho{}_s + \\rho{}_a|$")
ax[1].set_ylabel("$|v_s - v_a| / |v_s + v_a|$")
ax[1].set_xlabel("$r$ (AU)")
pl.tight_layout()
pl.savefig("figure_stable_solution_reldiff.png")
