import numpy as np
import matplotlib
#matplotlib.use("Agg")
import pylab as pl
import glob
import sys
import multiprocessing as mp
import os
import scipy.optimize as opt

pl.rcParams["figure.figsize"] = (12, 3)
pl.rcParams["text.usetex"] = True

gamma = 5. / 3.

au_in_si = 1.495978707e11 # m
yr_in_si = 365.25 * 24. * 3600. # s

def curve(x, A, B, C, D):
  return A + B * np.sin(C * x + D)

def jac_curve(x, A, B, C, D):
  J = np.zeros((len(x), 4))
  J[:, 0] = 1.
  J[:, 1] = np.sin(C * x + D)
  J[:, 2] = B * np.cos(C * x + D) * x
  J[:, 3] = B * np.cos(C * x + D)
  return J

def fit_curve(t, r, p0):
  fitt = []
  fitr = []
  for i in range(len(t)):
    if t[i] > 30. and t[i] < 60.:
      fitt.append(t[i])
      fitr.append(r[i])

  A,_ = opt.curve_fit(curve, fitt, fitr, p0 = p0,
                      jac = jac_curve)
  return A

def get_max_radius(f, num):
  ifile = open("lowres/" + f, 'r')
  timeline = ifile.readline()
  time = float(timeline.split()[2])
  ifile.close()
  data_lr = np.loadtxt("lowres/" + f)
  data_lr[:,4] = 1. - data_lr[:,4]
  data_lr[:,0] *= data_lr[:,4]
  data_lr_st = np.loadtxt("uhires/" + f)
  data_lr_st[:,4] = 1. - data_lr_st[:,4]
  data_lr_st[:,0] *= data_lr_st[:,4]
  data_hr = np.loadtxt("hires/" + f)
  data_hr[:,4] = 1. - data_hr[:,4]
  data_hr[:,0] *= data_hr[:,4]
  return num, time / yr_in_si, data_lr[:,0].max() / au_in_si, \
         data_hr[:,0].max() / au_in_si, data_lr_st[:,0].max() / au_in_si

pool = mp.Pool(16)
results = []
numsnap = 0
for f in sorted(glob.glob("hires/snapshot_*.txt")):
  num = int(f[-8:-4])
  filename = os.path.basename(f)
  if num > 0:
    results.append(pool.apply_async(get_max_radius, (filename, num,)))
  numsnap += 1

times = np.zeros(numsnap)
radius_low = np.zeros(numsnap)
radius_high = np.zeros(numsnap)
radius_st = np.zeros(numsnap)
for result in results:
  index, t, rlow, rhigh, rst = result.get()
  times[index] = t
  radius_low[index] = rlow
  radius_high[index] = rhigh
  radius_st[index] = rst
  print index, t, rlow, rhigh, rst

Alow = fit_curve(times, radius_low, (18.75, 0.5, 4.9, 0.))
Ahigh = fit_curve(times, radius_high, (18.75, 0.6, 4.6, 0.))

print "low:", 2.*np.pi / Alow[2]
print "high:", 2. * np.pi / Ahigh[2]

pl.plot(times, radius_low, label = "100 cells")
pl.plot(times, radius_high, label = "1000 cells")
pl.plot(times, radius_st, label = "2000 cells")
trange = np.arange(30., 60., 0.1)
pl.plot(trange, curve(trange, *Alow), ".")
pl.plot(trange, curve(trange, *Ahigh), ".")
pl.ylim(15., 25.)
pl.xlabel("time (yr)")
pl.ylabel("ionization radius (AU)")
pl.title(r"$\Delta{}t = 3.24 \times{} 10^{-4}$ yr")
pl.legend(loc = "best")
pl.tight_layout()
#pl.savefig("ionization_radius.png")
pl.show()
