#!/usr/bin/env python
from __future__ import print_function
from mpi4py import MPI
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
import time
import numpy as np
import sys
import os
import pandas as pd
import tools
import wrapper
# import scipy.stats

time.sleep(5)

t0 = time.time()

###########################################
#  Globals ##############################
###########################################

# PATH GLOBALS & tests to check they are as expected

# Save this path (expect data and other needed executables present in relative paths!!!)
THIS_PATH = os.path.dirname(__file__)

# Output path: it will be in scratch
# OUTPATH = os.path.abspath(os.path.join(os.path.expanduser('~'), 'scratch', 'output'))
OUTPATH = os.path.join(os.environ['GROUP_HOME'], "output")
outpk = os.path.join(os.environ['GROUP_HOME'], "failed_Pk")
outpktemp = os.path.join(outpk, "temp")

if not os.path.isdir(OUTPATH):
    raise Exception(OUTPATH + ' not there!')

# The CLASS Boltzmann code installation location
# Expects the CLASS code compiled here!!!
CLASSPATH = os.path.join(THIS_PATH, 'class')

# How to find the EFT model files and tools
# Expects the EFT code compiled here!!!
EFT_PATH = os.path.join(THIS_PATH, 'RedshiftBiasEFTwithFFT')

# Import the wrapper module. We assume it's in the same directory as here
# if os.path.isdir(EFT_PATH):
#     import sys
#     if EFT_PATH not in sys.path:
#         sys.path.append(EFT_PATH)
#     import wrapper
# else:
#     raise Exception('Module not found at ' + EFT_PATH)


# COSMOLOGICAL GLOBALS: fiducial model (should match input sim data!)
dfcosmo = pd.read_csv(os.path.join(THIS_PATH, 'input',
                                   'DataFrameCosmosims.csv'), index_col=0)
# allsims = ['LightConeDida', 'ChallengeA', 'ChallengeJapanCMASS2', 'LightConeHectorPatchyTrue']
simtype = 'ChallengeA'
series_cosmo = dfcosmo.loc[simtype]
# allgrids = ['ChallengeFullHDvFFT_IR06', 'ChallengeJapanHDvFFT_IR06',
#             'LightConeDidaHDvFFT_IR06', 'ChallengeHDvFFT_IR06',
#             'LightConeHectorPatchyz55HDvFFT_IR06']
gridname = 'GridSDSSChallenge'

# Cosmological parameters
ocfid = series_cosmo.loc['Omega_m'] * \
    series_cosmo.loc['h']**2 - series_cosmo.loc['omega_b']
obfid = series_cosmo.loc['omega_b']
fb = obfid/ocfid
Omega_mfid = series_cosmo['Omega_m']
h_TRUE = series_cosmo['h']
lnAs_TRUE = series_cosmo['lnAs']
z_pk = series_cosmo['z_pk']
nsfid = series_cosmo['ns']


###########################################
#  Grid ##############################
###########################################

lnAsmin, lnAsmax, nbinsAs = (2.5, 3.7, 120)
Ommin, Ommax, nbinsOm = (0.2, 0.42, 55)
hmin, hmax, nbinsh = (0.58, 0.78, 60)
lnAs = np.linspace(lnAsmin, lnAsmax, nbinsAs)
Om = np.linspace(Ommin, Ommax, nbinsOm)
h = np.linspace(hmin, hmax, nbinsh)
# The following gives me a matrix N x 3, with h varying fastest, then Om, then lnAs
thetatab = np.array(np.meshgrid(lnAs, Om, h, indexing='ij')).reshape(3, -1).T

###########################################
#  Functions  ###########################
###########################################

# The pre-computed functions f(\Omega_m) and \Omega_m(f) (capital \Omega, not \omega which depends on h)


# def cH(Om, a):
#     """Hubble in conformal time"""
#     return np.sqrt(Om/a+a**2*(1-Om))


# def DgN(Om, a):
#     """Linear growth factor"""
#     return 5./2*Om*cH(Om, a)/a*scipy.integrate.quad(lambda x: cH(Om, x)**-3, 0, a)[0]


def get_EFT_pars(lnAs, Om, h, config_file, nrank):
    omega_m = Om*h**2
    pars = tools.get_config(bigconfig_file=config_file, cat=True)
    pars['omega_cdm'] = omega_m / (1. + fb)
    pars['omega_b'] = omega_m * fb / (1. + fb)
    pars['h'] = h
    pars['ln10^{10}A_s'] = lnAs
    return pars

# def wR(x):
#     """Window function, FT of a top-hat"""
#     return (3*(-(x*np.cos(x)) + np.sin(x)))/x**3

# Da sistemare!!!


def CompPterms(theta, config_file, nrank):
    """Run the power spectra specified parameters.
    Assume PS are given in a file concatenated by multipoles.
    Inputs:
    - theta: Set of parameters for a single power spectrum
    - nrank: Rank of the MPI process
    - run: Index of the run
    Outputs:
    Plin, Ploop"""
    lnAs, Om, h = theta
    # Get the EFT parameters
    pars = get_EFT_pars(lnAs, Om, h, config_file, nrank)
    # Runs Pierre's code and save output to folder in path.
    path = wrapper.run_zbEFT(pars)
    # Get the k-values
    Plin = np.loadtxt(os.path.abspath(
        os.path.join(path, 'PowerSpectraLinear.dat')))
    Ploop = np.loadtxt(os.path.abspath(
        os.path.join(path, 'PowerSpectra1loop.dat')))
    # np.save(os.path.join(OUTPATH, 'textfiles/TablePloopSC%shalf%srun%sstep%s' %
    #                      (gridname, halfnum, 250*run+nrank, i)), Ploop)
    # np.save(os.path.join(OUTPATH, 'textfiles/TablePlinSC%shalf%srun%sstep%s' %
    #                      (gridname, halfnum, 250*run+nrank, i)), Plin)
    # # np.save(os.path.join(OUTPATH,'textfiles/Tables8SC%shalf%srun%sstep%s'%(gridname,halfnum,250*run+nrank,i)),s8)
    # # np.save(os.path.join(OUTPATH,'textfiles/TablepkClassSChallAhalf%srun%sstep%s'%(halfnum,250*run+nrank,i)),respkC)
    # np.save(os.path.join(OUTPATH, 'textfiles/TablecoordSC%shalf%srun%sstep%s' %
    #                      (gridname, halfnum, 250*run+nrank, i)), theta)
    return Plin, Ploop


###########################################
#  Data  ###########################
###########################################

totruns = 20
totranks = 25
Ntheta = len(thetatab)  # Total points of the full grid
lenrun = int(Ntheta / totruns)  # Length of chunks in the first run
sizerank = lenrun / totranks  # Length of sub-chunks processed by a single core

# Vector of failed points, (run, rank) tuples
allfailed = np.load("allfailed.npy")
# This are the total jobs we submit through the bash script
njobs = len(allfailed)
ncores = size  # This is now a small number of processors
# the sub-chunks are now subdivided in more cores
sizered = int(sizerank / ncores)

failedidx = int(sys.argv[1])  # Argument that tells us which index we are processing

nrunfail, nrankfail = allfailed[failedidx]
idxstart = nrunfail * lenrun + nrankfail * sizerank
idxend = nrunfail * lenrun + (nrankfail + 1) * sizerank
thetarun = thetatab[idxstart:idxend]  # Failed piece of the first run
# Small piece of the subchunk
smaller = thetarun[rank*sizered:(rank+1)*sizered]
allPlin = []
allPloop = []
for j, theta in enumerate(smaller):
    idx = idxstart + rank * sizered + j
    Plin, Ploop = CompPterms(theta, 'config.ini', rank)
    idxcol = np.full([Plin.shape[0], 1], idx)
    allPlin.append(np.hstack([Plin, idxcol]))
    allPloop.append(np.hstack([Ploop, idxcol]))
    if (j == 0) or ((j+1) % 50 == 0):
        print("Check parameters: ", thetatab[idx], theta)
        np.save(os.path.join(outpktemp, "Plin_run%s_rank%s_newrank%si%s.npy" %
                             (str(nrunfail), str(nrankfail), str(rank), str(j))), np.array(allPlin))
        np.save(os.path.join(outpktemp, "Ploop_run%s_rank%s_newrank%si%s.npy" %
                             (str(nrunfail), str(nrankfail), str(rank), str(j))), np.array(allPloop))
np.save(os.path.join(outpk, "Plin_run%s_rank%s_newrank%s.npy" %
                     (str(nrunfail), str(nrankfail), str(rank))), np.array(allPlin))
np.save(os.path.join(outpk, "Ploop_run%s_rank%s_newrank%s.npy" %
                     (str(nrunfail), str(nrankfail), str(rank))), np.array(allPloop))

print(time.time() - t0)
