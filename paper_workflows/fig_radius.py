import numpy as np
import matplotlib
matplotlib.use("Agg")
import pylab as pl
import glob
import sys

if len(sys.argv) < 2:
  print "Usage python fig_radius.py amplitude"
  exit()

amplitude = float(sys.argv[1])

pl.rcParams["text.usetex"] = True
pl.rcParams["figure.figsize"] = (6, 4)
pl.rcParams["font.size"] = 14
pl.rcParams["axes.labelsize"] = 18

au_in_si = 1.496e11
yr_in_si = (365.25 * 24. * 3600.)

widths = [5]

fig, ax = pl.subplots(1, 1, sharex = True, sharey = True)
ax = [ax]

for i in [0]:
  for res in [2700]:
    file = "convergence_seed_w{0}_{1}_{2}{3}_radius.dat".format(
      widths[i], res, 'p' if amplitude > 0. else 'm', abs(amplitude))
    print "Plotting", file, "..."

    fp = np.memmap(file, dtype = 'd', mode = 'r')
    data = fp.reshape((-1, 3))

    ax[i].plot(data[:,0] / yr_in_si, data[:,1] / au_in_si,
               label = "{0} cells".format(res))
  ax[i].axhline(y = 30, linestyle = "--", color = "k")
  ax[i].axhline(y = 10, linestyle = "-", color = "k")
  ax[i].set_title("$W = {0}$ AU".format(widths[i]))
  ax[i].set_ylabel("$R_I$ (AU)")

ax[-1].set_xlabel("$t$ (yr)")
pl.tight_layout()
pl.savefig("fig_radius_{sign}{amplitude}.png".format(
  sign = 'p' if amplitude > 0. else 'm', amplitude = abs(amplitude)))
