import numpy as np
import matplotlib
matplotlib.use("Agg")
import pylab as pl
import scipy.special.lambertw as lambertw
import mpl_toolkits.axes_grid1.inset_locator as inset

pl.rcParams["text.usetex"] = True
pl.rcParams["figure.figsize"] = (4.5, 5)

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

def plot(f, ax, zoomrho, zoomv):
  data = np.loadtxt("convergence_stable_" + f)
  data[:,0] /= au_in_si
  data[:,1] *= 0.001
  data[:,2] *= 0.001
  W = float(f[1:2])
  label = None
  if W < 3.:
    label = "$W = {0:.0f}$ AU".format(W)
  ax[0].semilogy(data[:,0], data[:,1], "-", label = label)
  label = None
  if W > 2.:
    label = "$W = {0:.0f}$ AU".format(W)
  ax[1].plot(data[:,0], data[:,2], "-", label = label)
  zoomrho.semilogy(data[:,0], data[:,1], "-")
  zoomv.plot(data[:,0], data[:,2], "-")

fig, ax = pl.subplots(2, 1, sharex = "col")
zoomrho = inset.inset_axes(ax[0], width = "40%", height = "50%", loc = 1)
zoomv = inset.inset_axes(ax[1], width = "100%", height = "100%", loc = 4,
                         bbox_to_anchor = (0.6, 0.1, 0.4, 0.5),
                         bbox_transform = ax[1].transAxes)

plot("w1_2700.txt", ax, zoomrho, zoomv)
plot("w2_2700.txt", ax, zoomrho, zoomv)
plot("w3_2700.txt", ax, zoomrho, zoomv)
plot("w4_2700.txt", ax, zoomrho, zoomv)
plot("w5_2700.txt", ax, zoomrho, zoomv)

ax[0].plot(ra, rhoa, "k--", linewidth = 0.8)
ax[1].plot(ra, va, "k--", linewidth = 0.8)

zoomrho.plot(ra, rhoa, "k--", linewidth = 0.8)
zoomv.plot(ra, va, "k--", linewidth = 0.8)

ax[0].set_ylabel("$\\rho{}$ (g cm$^{-3}$)")
ax[1].set_ylabel("$v$ (km s$^{-1}$)")

ax[0].set_ylim(3.e-18, 2.e-16)
ax[1].set_ylim(-50., -8.)
zoomrho.set_xlim(25., 35.)
zoomv.set_xlim(25., 35.)
zoomrho.set_ylim(2.e-17, 4.2e-17)
zoomv.set_ylim(-33., -26.)
ax[1].set_xlabel("$r$ (AU)")

ax[0].legend(loc = "lower left", ncol = 2)
ax[1].legend(loc = "upper left", ncol = 2)

pl.tight_layout()
pl.savefig("fig_convergence_stable.eps", dpi = 300)
pl.close()
