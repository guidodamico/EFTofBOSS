import os
import numpy as np

nruns = 20
nranks = 25
ntot = nruns * nranks
pathpk = os.path.join(os.environ["SCRATCH"], "Pk")

linfailed = []
loopfailed = []
for i in range(nruns):
    for j in range(nranks):
        checklin = os.path.isfile(os.path.join(
            pathpk, "Plin_run%d_rank%d.npy" % (i, j)))
        checkloop = os.path.isfile(os.path.join(
            pathpk, "Ploop_run%d_rank%d.npy" % (i, j)))
        if not checklin:
            print("Failed linear run %d, rank %d" % (i, j))
            linfailed.append((i, j))
        if not checkloop:
            print("Failed loop run %d, rank %d" % (i, j))
            loopfailed.append((i, j))

print("Linear failed: %d over %d, %f %%" %
      (len(linfailed), ntot, 100 * float(len(linfailed)) / ntot))
print("Loop failed: %d over %d, %f %%" %
      (len(loopfailed), ntot, 100 * float(len(loopfailed)) / ntot))

np.save("linfailed.npy", np.array(linfailed))
np.save("loopfailed.npy", np.array(loopfailed))
