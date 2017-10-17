import numpy as np
import pylab as pl
import glob
import sys
import multiprocessing as mp

gamma = 5. / 3.

au_in_si = 1.495978707e11 # m
yr_in_si = 365.25 * 24. * 3600. # s

last = 10000
if len(sys.argv) > 1:
  last = int(sys.argv[1])

def plot(f):
  ifile = open(f, 'r')
  timeline = ifile.readline()
  time = float(timeline.split()[2])
  ifile.close()
  data = np.loadtxt(f)
#  refdata = np.loadtxt("../hires/" + f)
  fig, ax = pl.subplots(2, 2, sharex = "col")
  data[:,0] /= au_in_si
#  ax[0][0].plot(refdata[:,0], refdata[:,1], "r-")
#  ax[0][1].plot(refdata[:,0], refdata[:,2], "r-")
#  ax[1][0].plot(refdata[:,0], refdata[:,3], "r-")
  ax[0][0].plot(data[:,0], data[:,1]*0.001, "k-")
  ax[0][1].plot(data[:,0], data[:,2]*0.001, "k-")
  ax[1][0].plot(data[:,0], data[:,3], "k-")
  ax[1][1].plot(data[:,0], data[:,4], "k-")
  ax[0][0].set_ylim(0., 3.2e-10)
  ax[1][0].set_ylim(0., 2.2e-12)
  ax[1][1].set_ylim(-0.1, 1.1)
  ax[0][0].set_ylabel("density (g cm$^{-3}$)")
  ax[0][1].set_ylabel("velocity (km s$^{-1}$)")
  ax[1][0].set_ylabel("pressure (kg m$^{-1}$ s$^{-2}$)")
  ax[1][1].set_ylabel("neutral fraction")
  ax[1][0].set_xlabel("radius (AU)")
  ax[1][1].set_xlabel("radius (AU)")
  ax[0][1].set_title("t = {t:.2e} yr".format(t = time / yr_in_si))
  pl.tight_layout()
  pl.savefig("{name}.png".format(name = f[:-4]))
  pl.close()

  return f, data[:,1].max(), data[:,3].max()

pool = mp.Pool(16)
results = []
for f in sorted(glob.glob("snapshot_*.txt")):
  num = int(f[-7:-4])
  if num < last:
    results.append(pool.apply_async(plot, (f,)))

maxrho = 0.
maxP = 0.
for result in results:
  f, rho, P = result.get()
  maxrho = max(maxrho, rho)
  maxP = max(maxP, P)
  print f

print "maxrho:", maxrho
print "maxP:", maxP
