import os
import numpy as np

pathpk = os.path.join(os.environ["GROUP_HOME"], "Pk")
pathgrid = os.path.join(os.environ["GROUP_HOME"], "grids")
nruns = 20
nranks = 25
lenbatch = 792

gridlin = []
gridloop = []
for i in range(nruns):
    print("Run ", i)
    for j in range(nranks):
        Plin = np.load(os.path.join(
            pathpk, "Plin_run%s_rank%s.npy" % (str(i), str(j))))
        Ploop = np.load(os.path.join(
            pathpk, "Ploop_run%s_rank%s.npy" % (str(i), str(j))))
        gridlin.append(Plin[:, :, :-1])
        gridloop.append(Ploop[:, :, :-1])
        checklin = ((i * nranks + j) * lenbatch == int(Plin[0, 0, -1]))
        checkloop = ((i * nranks + j) * lenbatch == int(Ploop[0, 0, -1]))
        if not checklin:
            print(i, j, (i * nranks + j) * lenbatch, Plin[0, 0, -1])
        if not checkloop:
            print(i, j, (i * nranks + j) * lenbatch, Ploop[0, 0, -1])
p1=np.concatenate(gridlin)
np.save(os.path.join(pathgrid, "TablePlin_z0p5617_bins120x55x60.npy"), p1)
p2=np.concatenate(gridloop)
np.save(os.path.join(pathgrid, "TablePloop_z0p5617_bins120x55x60.npy"), p2)

print("I'm done!")
