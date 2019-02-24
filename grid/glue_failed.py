import os
import numpy as np

allfailed = np.load("allfailed.npy")
pathpk = os.path.join(os.environ["GROUP_HOME"], "failed_Pk")
pathnewpk = os.path.join(os.environ["GROUP_HOME"], "new_Pk")
njobs = 33

for i, tup in enumerate(allfailed):
    nrunfail, nrankfail = tup
    Plins = [os.path.join(pathpk, "Plin_run%s_rank%s_newrank%s.npy" %
                     (str(nrunfail), str(nrankfail), str(j))) for j in range(njobs)]
    Ploops = [os.path.join(pathpk, "Ploop_run%s_rank%s_newrank%s.npy" %
                     (str(nrunfail), str(nrankfail), str(j))) for j in range(njobs)]
    checklin = [os.path.isfile(f) for f in Plins]
    checkloop = [os.path.isfile(f) for f in Ploops]
    if (sum(checklin) == njobs) and (sum(checkloop) == njobs):
        print(i, nrunfail, nrankfail)
        fileslin = [np.load(p) for p in Plins]
        filesloop = [np.load(p) for p in Ploops]
        p1 = np.concatenate(fileslin)
        p2 = np.concatenate(filesloop)
        np.save(os.path.join(pathnewpk, "Plin_run%s_rank%s.npy" % (str(nrunfail), str(nrankfail))), p1)
        np.save(os.path.join(pathnewpk, "Ploop_run%s_rank%s.npy" % (str(nrunfail), str(nrankfail))), p2)

print("I'm done!")
