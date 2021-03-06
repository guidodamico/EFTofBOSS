#!/usr/bin/env python

###########################################
###  Imports  #############################
###########################################

from __future__ import print_function

import os
import mpi4py
mpi4py.rc.initialize = False
mpi4py.rc.finalize = False
#os.environ["TMPDIR"] = "/tmp"
import emcee
import numpy as np
from numpy import unravel_index
import scipy.stats
import scipy.optimize as op
from scipy.interpolate import interp1d
import pandas as pd
import os.path as opa
import time
from scipy import stats
import sys
print("imported")
sys.stdout.flush()
###########################################
###  Globals ##############################
###########################################

# PATH GLOBALS & tests to check they are as expected

# Save this path (expect data and other needed executables present in relative paths!!!)
THIS_PATH  =  opa.dirname(__file__)

# Import local packages (to be saved in same folder)
import APpowerspectraNkmu
import WindowFunctionFourier

# Data paths
INPATH = opa.abspath(opa.join(THIS_PATH,'input'))
#INPATH2 = opa.abspath('/scratch/users/kokron/')
OUTPATH = opa.abspath(opa.join(THIS_PATH,'output_bestfit')) 
#print(OUTPATH)
if not opa.isdir(OUTPATH): raise Exception(OUTPATH + ' not there!')


###########################################
###  Functions  ###########################
###########################################

### rescale number density
#nd = 1
#km = 1
#knl = 1
nd = 3e-4
km = 0.7
knl = 0.6
#nd = 20
#km = 0.7
#knl = 0.6
#####################################################################################################################################################
################## Sound horizon at decoupling #######################
###
### Option 'withPlanck' added (l.433) in lnlike: setup it up at l.518.
###
################################################### - Pierre, 10/08/18
#####################################################################################################################################################
# speed of light [km/s]
C = 299792.458 
# omega_gamma = Omega_gamma h^2: normalized physical photon density today (T_cmb = 2.7255 (CLASS))
OG = 2.47282e-5
#Nur: Number of ultra-relativistic species (CLASS):
NUR = 3.046 
# omega_radiation
ORAD = (1.+ NUR*7./8.*(4./11.)**(4./3.))*OG

# Baryon-photon decoupling redshift (PLANCK 2015 TT,TE,EE+lowP+lensing (Table 4)):
ZD = 1059.94
SIGMA_ZD = 0.30 # rd(zd+sigma)-dr(zd-sigma) < 0.2 sigma_rd: we take zd to be a delta function
# Sound horizon at decoupling [Mpc] (PLANCK 2015 TT,TE,EE+lowP+lensing (Table 4)):
RD = 147.09
SIGMA_RD = 0.26

# Baryon-photon decoupling redshift (PLANCK 2015 TT,TE,EE+lowP+lensing (Table 4)):
#ZD = 1059.62
#SIGMA_ZD = 0.31 # rd(zd+sigma)-dr(zd-sigma) < 0.2 sigma_rd: we take zd to be a delta function
# Sound horizon at decoupling [Mpc] (PLANCK 2015 TT,TE,EE+lowP+lensing (Table 4)):
#RD = 147.41
#SIGMA_RD = 0.30

def rs(Om,h,f_fid):
    om = Om*h**2
    ob = om * f_fid/(f_fid+1.) # f_fid: fiducial ratio omega_b/omega_c
    R = 0.75 * ob/OG
    result = 2.*C/100./ np.sqrt(3.*R*om)* np.log( ( np.sqrt(1.+ZD+R) + np.sqrt((1.+ZD)*R*ORAD/om+R) ) / np.sqrt(1.+ZD) / (1.+np.sqrt(R*ORAD/om)) )
    return result
#####################################################################################################################################################


def embed_Pi(ploop, masktriangle):
    #Embed the Pi into a bigger shape that has appropriate prepadding and postpadding so the bisp term can be easily added on
    #Dimension of Pi is [nterms, 3, 100]
    #return with dimension [nterms+1, 4, 100+nkbisp]
    nkbisp = sum(masktriangle)
    nkp = ploop.shape[1]
    #Cast into bigger array of dim p
    big_array = np.zeros(shape=(ploop.shape[0]+1,100+nkbisp)) 
    return big_array
def get_Pi_for_marg(Ploop,kfull,b1,bisp=None):

    nk = len(kfull)
    Onel0 = np.array([np.ones(nk),np.zeros(nk),np.zeros(nk)])
    kl0 = np.array([kfull,np.zeros(nk),np.zeros(nk)])
    kl2 = np.array([np.zeros(nk),kfull,np.zeros(nk)])

    Pi = np.array([Ploop[:,3,:]+b1*Ploop[:,7,:],
             (Ploop[:,15,:]+b1*Ploop[:,12,:]) / knl**2,
             (Ploop[:,16,:]+b1*Ploop[:,13,:]) / km**2,
             (Ploop[:,17,:]+b1*Ploop[:,14,:]) / km**2,
             Onel0 / nd,
             # PCA b9 b10
             (kl0**2 + kl2**2) / nd / km**2,
             kl2**2 / nd / km**2])

    if withBisp:
        #b8 is not marginalized with bisp but b11 is.
        #print(Ploop[:,3,:].shape, Ploop[:,7,:].shape, postpad.shape)

        #don't have b8 term

        Pi = np.array([Ploop[:,3,:]+b1*Ploop[:,7,:],
                 (Ploop[:,15,:]+b1*Ploop[:,12,:]) / knl**2,
                 (Ploop[:,16,:]+b1*Ploop[:,13,:]) / km**2,
                 (Ploop[:,17,:]+b1*Ploop[:,14,:]) / km**2,
                 (kl0**2 + kl2**2) / nd / km**2,
                 kl2**2 / nd / km**2])
    return Pi
    
    

def get_Covbi_for_marg(Pi_data,Cinv,sigma=200):

    Covbi = np.dot(Pi_data,np.dot(Cinv,Pi_data.T))+ 1./sigma**2*np.identity(Pi_data.shape[0])
    return Covbi
    

def Hubble(Om,z):
    return ((Om)*(1+z)**3.+(1-Om))**0.5

def DA(Om,z):
    r = scipy.integrate.quad(lambda x:1./Hubble(Om,x), 0, z)[0]
    return r/(1+z)  

def get_AP_param(Om,Om_fid):
        
    qperp  =  DA(Om,z_pk)/DA(Om_fid,z_pk)
    qpar  =  Hubble(Om_fid,z_pk)/Hubble(Om,z_pk)
    
    return qperp,qpar

def check_if_multipoles_k_array(setk):
    return setk[len(setk)/3] == setk[0]    
            


def get_grid(gridname,nbinsAs=100,nbinsOm =48,nbinsh=48,withBisp=False):
    
    """ Computes the power spectra given the b_i and the EFT power spectra
        Inputs
        ------
        gridname : The name of grid associated to the sim
        nbinsAs : number of bins for As (default is 100)
        nbinsAs : number of bins for Om and h (default is 50)
        withBisp : whether or not to load grid for bispectrum (only works for simtype = LightConeHector)
        Outputs
        ------
        The min,max values for the three parameters as well as the interpolation for the linear and loop power spectra
    """
    
    thetatab = np.load(opa.abspath(opa.join(INPATH,'GridsEFT/Tablecoord%s.npy'%gridname)))

    theta3D = thetatab.reshape((nbinsAs,nbinsOm,nbinsh,3))

    lnAstab = theta3D[:,0,0,0]
    Omtab = theta3D[0,:,0,1]
    htab = theta3D[0,0,:,2]


    lnAsmin = lnAstab.min()
    lnAsmax = lnAstab.max()
    Ommin = Omtab.min()
    Ommax = Omtab.max()
    hmin = htab.min()
    hmax = htab.max()

    TablePlin = np.load(opa.abspath(opa.join(INPATH,'GridsEFT/TablePlin%s.npy'%gridname)))
    TablePloop = np.load(opa.abspath(opa.join(INPATH,'GridsEFT/TablePloop%s.npy'%gridname)))
    #Tablesigsq = np.load(opa.abspath(opa.join(INPATH,'GridsEFT/Tablesigsq%s.npy'%gridname)))

    Plininterp = scipy.interpolate.RegularGridInterpolator((lnAstab,Omtab,htab),TablePlin.reshape((nbinsAs,nbinsOm,nbinsh,TablePlin.shape[-2],TablePlin.shape[-1])))
    Ploopinterp = scipy.interpolate.RegularGridInterpolator((lnAstab,Omtab,htab),TablePloop.reshape((nbinsAs,nbinsOm,nbinsh,TablePloop.shape[-2],TablePloop.shape[-1])))
    #Sigsqinterp = scipy.interpolate.RegularGridInterpolator((lnAstab,Omtab,htab),Tablesigsq.reshape((nbinsAs,nbins,nbins)))
    
    #interpolations = [Plininterp,Ploopinterp,Sigsqinterp]
    interpolations = [Plininterp,Ploopinterp]
    if withBisp:
        TableBisp = np.load(opa.abspath(opa.join(INPATH,'GridsEFT/TableBisp%s.npy'%gridname)))
        Bispinterp = scipy.interpolate.RegularGridInterpolator((lnAstab,Omtab,htab),TableBisp.reshape((nbinsAs,nbinsOm,nbinsh,TableBisp.shape[-2],TableBisp.shape[-1])))
        #interpolations = [Plininterp,Ploopinterp,Sigsqinterp,Bispinterp]
        interpolations = [Plininterp,Ploopinterp,Bispinterp]

        
    return lnAsmin,lnAsmax,Ommin,Ommax,hmin,hmax,interpolations

################################

################################

################################

def computePS(cvals,datalin,dataloop,setkin,setkout,sigsq=0):
    
    """ Computes the power spectra given the b_i and the EFT power spectra
        Inputs
        ------
        bvals : The values for the b_i
        datalin : the linear power spectra from the EFT, with shape (multipoles, b_i, k)
        dataloop : the loop power spectra from the EFT, with shape (multipoles, b_i, k)
        setkin : the values of k from the input power spectra (must match datalin and dataloop)
        setkout : the values of k for the output power spectra
        sigsq: The sigma^2 that needs to be removed, corresponding to the constant piece of the loop terms
        Outputs
        ------
        The power spectra multipoles, non-concatenated
    """
    
    datalin0,datalin2,datalin4 = datalin
    data0,data2,data4 = dataloop
    b1,c2,b3,c4,b5,b6,b7,b8,b9,b10 = cvals#

    # Add a PCA between b2 and b4 to disantangle the degeneracy:
    b2 = 0.5 * (c2 + c4)
    b4 = 0.5 * (c2 - c4)

    
    # the columns of the Ploop data files.
    cvals = np.array([1,b1,b2,b3,b4,b1*b1,b1*b2,b1*b3,b1*b4,b2*b2,b2*b4,b4*b4, b1*b5/knl**2, b1*b6/km**2, b1*b7/km**2, b5/knl**2, b6/km**2 ,b7/km**2, b8*nd/km**2 ,b9*nd/km**2, b10*nd/km**2])
    #cvals = np.array([1./(2.*np.pi),b1/(2.*np.pi),b2/(2.*np.pi),b3/(2.*np.pi),b4/(2.*np.pi),b1*b1/(2.*np.pi),b1*b2/(2.*np.pi),b1*b3/(2.*np.pi),b1*b4/(2.*np.pi),b2*b2/(2.*np.pi),b2*b4/(2.*np.pi),b4*b4*2/(4.+np.pi),b1*b5,b1*b6,b1*b7,b5,b6 ,b7,b8 ,b9,b10])
    
    # Check the k-arrays are in the right format (not concatenated for multipoles)
    if check_if_multipoles_k_array(setkin):
        setkin = setkin[:len(setkin)/3]
    if check_if_multipoles_k_array(setkout):
        setkout = setkin[:len(setkout)/3]
         
        
    P0 = interp1d(setkin,np.dot(cvals,data0)+datalin0[0]+b1*datalin0[1]+b1*b1*datalin0[2] - 2*(-b1 + b2 + b4)**2*sigsq)(setkout)
    P2 = interp1d(setkin,np.dot(cvals,data2)+datalin2[0]+b1*datalin2[1]+b1*b1*datalin2[2])(setkout)
    P4 = interp1d(setkin,np.dot(cvals,data4)+datalin4[0]+b1*datalin4[1]+b1*b1*datalin4[2])(setkout)
    
    return np.array([P0,P2,P4])

def computeLoop(cvals,datalin,dataloop,setkin,setkout,sigsq=0):
    
    datalin0,datalin2,datalin4 = datalin
    data0,data2,data4 = dataloop
    b1,c2,b3,c4,b5,b6,b7,b8,b9,b10 = cvals

    # Add a PCA between b2 and b4 to disantangle the degeneracy:
    b2 = 0.5 * (c2 + c4)
    b4 = 0.5 * (c2 - c4)

    # the columns of the Ploop data files.
    cvals = np.array([1,b1,b2,b3,b4,b1*b1,b1*b2,b1*b3,b1*b4,b2*b2,b2*b4,b4*b4, b1*b5/knl**2, b1*b6/km**2, b1*b7/km**2, b5/knl**2, b6/km**2 ,b7/km**2, b8*nd/km**2 ,b9*nd/km**2, b10*nd/km**2])
    #cvals = np.array([1./(2.*np.pi),b1/(2.*np.pi),b2/(2.*np.pi),b3/(2.*np.pi),b4/(2.*np.pi),b1*b1/(2.*np.pi),b1*b2/(2.*np.pi),b1*b3/(2.*np.pi),b1*b4/(2.*np.pi),b2*b2/(2.*np.pi),b2*b4/(2.*np.pi),b4*b4*2/(4.+np.pi),b1*b5,b1*b6,b1*b7,b5,b6 ,b7,b8 ,b9,b10])
    
    # Check the k-arrays are in the right format (not concatenated for multipoles)
    if check_if_multipoles_k_array(setkin):
        setkin = setkin[:len(setkin)/3]
    if check_if_multipoles_k_array(setkout):
        setkout = setkin[:len(setkout)/3]
         
        
    P0 = interp1d(setkin,np.dot(cvals,data0))(setkout)
    P2 = interp1d(setkin,np.dot(cvals,data2))(setkout)
    P4 = interp1d(setkin,np.dot(cvals,data4))(setkout)


def computeLoop2(cvals,datalin,dataloop,setkin,setkout,sigsq=0):
    
    datalin0,datalin2,datalin4 = datalin
    data0,data2,data4 = dataloop
    b1,c2,b3,c4,b5,b6,b7,b8,b9,b10 = cvals#

    # Add a PCA between b2 and b4 to disantangle the degeneracy:
    b2 = 0.5 * (c2 + c4)
    b4 = 0.5 * (c2 - c4)

    
    # the columns of the Ploop data files.
    cvals = np.array([0,b1,b2,b3,b4,b1*b1,b1*b2,b1*b3,b1*b4,b2*b2,b2*b4,b4*b4, b1*b5/knl**2, b1*b6/km**2, b1*b7/km**2, b5/knl**2, b6/km**2 ,b7/km**2, b8*nd/km**2 ,b9*nd/km**2, b10*nd/km**2])
    #cvals = np.array([1./(2.*np.pi),b1/(2.*np.pi),b2/(2.*np.pi),b3/(2.*np.pi),b4/(2.*np.pi),b1*b1/(2.*np.pi),b1*b2/(2.*np.pi),b1*b3/(2.*np.pi),b1*b4/(2.*np.pi),b2*b2/(2.*np.pi),b2*b4/(2.*np.pi),b4*b4*2/(4.+np.pi),b1*b5,b1*b6,b1*b7,b5,b6 ,b7,b8 ,b9,b10])
    
    # Check the k-arrays are in the right format (not concatenated for multipoles)
    if check_if_multipoles_k_array(setkin):
        setkin = setkin[:len(setkin)/3]
    if check_if_multipoles_k_array(setkout):
        setkout = setkin[:len(setkout)/3]
         
        
    P0 = interp1d(setkin,np.dot(cvals,data0)+datalin0[0]+b1*datalin0[1]+b1*b1*datalin0[2] - 2*(-b1 + b2 + b4)**2*sigsq)(setkout)
    P2 = interp1d(setkin,np.dot(cvals,data2)+datalin2[0]+b1*datalin2[1]+b1*b1*datalin2[2])(setkout)
    P4 = interp1d(setkin,np.dot(cvals,data4)+datalin4[0]+b1*datalin4[1]+b1*b1*datalin4[2])(setkout)
    
    return np.array([P0,P2,P4])


def computeStoch(cvals,datalin,dataloop,setkin,setkout,sigsq=0):
    
    datalin0,datalin2,datalin4 = datalin
    data0,data2,data4 = dataloop
    b1,c2,b3,c4,b5,b6,b7,b8,b9,b10 = cvals

    # Add a PCA between b2 and b4 to disantangle the degeneracy:
    b2 = 0.5 * (c2 + c4)
    b4 = 0.5 * (c2 - c4)

    # the columns of the Ploop data files.
    cvals = np.array([0,b1,b2,b3,b4,b1*b1,b1*b2,b1*b3,b1*b4,b2*b2,b2*b4,b4*b4, b1*b5/knl**2, b1*b6/km**2, b1*b7/km**2, b5/knl**2, b6/km**2 ,b7/km**2, b8*nd/km**2 ,b9*nd/km**2, b10*nd/km**2])
    #cvals = np.array([1./(2.*np.pi),b1/(2.*np.pi),b2/(2.*np.pi),b3/(2.*np.pi),b4/(2.*np.pi),b1*b1/(2.*np.pi),b1*b2/(2.*np.pi),b1*b3/(2.*np.pi),b1*b4/(2.*np.pi),b2*b2/(2.*np.pi),b2*b4/(2.*np.pi),b4*b4*2/(4.+np.pi),b1*b5,b1*b6,b1*b7,b5,b6 ,b7,b8 ,b9,b10])
    
    # Check the k-arrays are in the right format (not concatenated for multipoles)
    if check_if_multipoles_k_array(setkin):
        setkin = setkin[:len(setkin)/3]
    if check_if_multipoles_k_array(setkout):
        setkout = setkin[:len(setkout)/3]
         
        
    P0 = interp1d(setkin,np.dot(cvals,data0)+datalin0[0]+b1*datalin0[1]+b1*b1*datalin0[2] - 2*(-b1 + b2 + b4)**2*sigsq)(setkout)
    P2 = interp1d(setkin,np.dot(cvals,data2)+datalin2[0]+b1*datalin2[1]+b1*b1*datalin2[2])(setkout)
    P4 = interp1d(setkin,np.dot(cvals,data4)+datalin4[0]+b1*datalin4[1]+b1*b1*datalin4[2])(setkout)

################################
################################
################################



'''
def computePS(cvals,datalin,dataloop,setkin,setkout,sigsq=0):
    
    """ Computes the power spectra given the b_i and the EFT power spectra
        Inputs
        ------
        bvals : The values for the b_i
        datalin : the linear power spectra from the EFT, with shape (multipoles, b_i, k)
        dataloop : the loop power spectra from the EFT, with shape (multipoles, b_i, k)
        setkin : the values of k from the input power spectra (must match datalin and dataloop)
        setkout : the values of k for the output power spectra
        sigsq: The sigma^2 that needs to be removed, corresponding to the constant piece of the loop terms
        Outputs
        ------
        The power spectra multipoles, non-concatenated
    """
    
    datalin0,datalin2,datalin4 = datalin
    data0,data2,data4 = dataloop[:,:18,:]
    b1,c2,b3,c4,b5,b6,b7,b8,b9,b10 = cvals

    # Add a PCA between b2 and b4 to disantangle the degeneracy:
    b2 = 0.5 * (c2 + c4)
    b4 = 0.5 * (c2 - c4)

    
    # the columns of the Ploop data files.
    cvals = np.array([1,b1,b2,b3,b4,b1*b1,b1*b2,b1*b3,b1*b4,b2*b2,b2*b4,b4*b4, b1*b5/knl**2, b1*b6/km**2, b1*b7/km**2, b5/knl**2, b6/km**2 ,b7/km**2])

    # Check the k-arrays are in the right format (not concatenated for multipoles)
    if check_if_multipoles_k_array(setkin):
        setkin = setkin[:len(setkin)/3]
    if check_if_multipoles_k_array(setkout):
        setkout = setkin[:len(setkout)/3]
         
        
    P0 = interp1d(setkin, np.dot(cvals,data0) + datalin0[0]+b1*datalin0[1]+b1*b1*datalin0[2])(setkout) + b8/nd + b9/nd/km**2 * setkout**2
    P2 = interp1d(setkin, np.dot(cvals,data2) + datalin2[0]+b1*datalin2[1]+b1*b1*datalin2[2])(setkout) + b9/nd/km**2 * setkout**2 + b10/nd/km**2 * setkout**2
    P4 = interp1d(setkin, np.dot(cvals,data4) + datalin4[0]+b1*datalin4[1]+b1*b1*datalin4[2])(setkout)
    
    return np.array([P0,P2,P4])
'''


def computeLin(cvals,datalin,dataloop,setkin,setkout,sigsq=0):
    
    """ Computes the power spectra given the b_i and the EFT power spectra
        Inputs
        ------
        bvals : The values for the b_i
        datalin : the linear power spectra from the EFT, with shape (multipoles, b_i, k)
        dataloop : the loop power spectra from the EFT, with shape (multipoles, b_i, k)
        setkin : the values of k from the input power spectra (must match datalin and dataloop)
        setkout : the values of k for the output power spectra
        sigsq: The sigma^2 that needs to be removed, corresponding to the constant piece of the loop terms
        Outputs
        ------
        The power spectra multipoles, non-concatenated
    """
    
    datalin0,datalin2,datalin4 = datalin
    data0,data2,data4 = dataloop
    b1,c2,b3,c4,b5,b6,b7,b8,b9,b10 = cvals

    # Add a PCA between b2 and b4 to disantangle the degeneracy:
    b2 = 0.5 * (c2 + c4)
    b4 = 0.5 * (c2 - c4)


    # Check the k-arrays are in the right format (not concatenated for multipoles)
    if check_if_multipoles_k_array(setkin):
        setkin = setkin[:len(setkin)/3]
    if check_if_multipoles_k_array(setkout):
        setkout = setkin[:len(setkout)/3]
         
        
    P0 = interp1d(setkin, datalin0[0]+b1*datalin0[1]+b1*b1*datalin0[2])(setkout)
    P2 = interp1d(setkin, datalin2[0]+b1*datalin2[1]+b1*b1*datalin2[2])(setkout)
    P4 = interp1d(setkin, datalin4[0]+b1*datalin4[1]+b1*b1*datalin4[2])(setkout)
    
    return np.array([P0,P2,P4])
'''
def computeLoop(cvals,datalin,dataloop,setkin,setkout,sigsq=0):
    
    datalin0,datalin2,datalin4 = datalin
    data0,data2,data4 = dataloop[:,:18,:]
    b1,c2,b3,c4,b5,b6,b7,b8,b9,b10 = cvals

    # Add a PCA between b2 and b4 to disantangle the degeneracy:
    b2 = 0.5 * (c2 + c4)
    b4 = 0.5 * (c2 - c4)

    cvals = np.array([1,b1,b2,b3,b4,b1*b1,b1*b2,b1*b3,b1*b4,b2*b2,b2*b4,b4*b4, b1*b5/knl**2, b1*b6/km**2, b1*b7/km**2, b5/knl**2, b6/km**2 ,b7/km**2])


    # Check the k-arrays are in the right format (not concatenated for multipoles)
    if check_if_multipoles_k_array(setkin):
        setkin = setkin[:len(setkin)/3]
    if check_if_multipoles_k_array(setkout):
        setkout = setkin[:len(setkout)/3]
         
        
    P0 = interp1d(setkin, np.dot(cvals,data0))(setkout) + b8/nd + b9/nd/km**2 * setkout**2
    P2 = interp1d(setkin, np.dot(cvals,data2))(setkout) + b9/nd/km**2 * setkout**2 + b10/nd/km**2 * setkout**2
    P4 = interp1d(setkin, np.dot(cvals,data4))(setkout)
    
    return np.array([P0,P2,P4])


def computeLoop2(cvals,datalin,dataloop,setkin,setkout,sigsq=0):
    
    datalin0,datalin2,datalin4 = datalin
    data0,data2,data4 = dataloop[:,:18,:]
    b1,c2,b3,c4,b5,b6,b7,b8,b9,b10 = cvals

    # Add a PCA between b2 and b4 to disantangle the degeneracy:
    b2 = 0.5 * (c2 + c4)
    b4 = 0.5 * (c2 - c4)

    cvals = np.array([0,b1,b2,b3,b4,b1*b1,b1*b2,b1*b3,b1*b4,b2*b2,b2*b4,b4*b4, b1*b5/knl**2, b1*b6/km**2, b1*b7/km**2, b5/knl**2, b6/km**2 ,b7/km**2])


    # Check the k-arrays are in the right format (not concatenated for multipoles)
    if check_if_multipoles_k_array(setkin):
        setkin = setkin[:len(setkin)/3]
    if check_if_multipoles_k_array(setkout):
        setkout = setkin[:len(setkout)/3]
         
        
    P0 = interp1d(setkin, np.dot(cvals,data0))(setkout) + b8/nd + b9/nd/km**2 * setkout**2
    P2 = interp1d(setkin, np.dot(cvals,data2))(setkout) + b9/nd/km**2 * setkout**2 + b10/nd/km**2 * setkout**2
    P4 = interp1d(setkin, np.dot(cvals,data4))(setkout)
    
    return np.array([P0,P2,P4])

def computeStoch(cvals,datalin,dataloop,setkin,setkout,sigsq=0):
    
    datalin0,datalin2,datalin4 = datalin
    data0,data2,data4 = dataloop[:,:18,:]
    b1,c2,b3,c4,b5,b6,b7,b8,b9,b10 = cvals

    # Add a PCA between b2 and b4 to disantangle the degeneracy:
    b2 = 0.5 * (c2 + c4)
    b4 = 0.5 * (c2 - c4)

    cvals = np.array([0,b1,b2,b3,b4,b1*b1,b1*b2,b1*b3,b1*b4,b2*b2,b2*b4,b4*b4, b1*b5/knl**2, b1*b6/km**2, b1*b7/km**2, b5/knl**2, b6/km**2 ,b7/km**2])


    # Check the k-arrays are in the right format (not concatenated for multipoles)
    if check_if_multipoles_k_array(setkin):
        setkin = setkin[:len(setkin)/3]
    if check_if_multipoles_k_array(setkout):
        setkout = setkin[:len(setkout)/3]
         
        
    P0 = interp1d(setkin, np.dot(cvals,data0))(setkout) + b8/nd + b9/nd/km**2 * setkout**2
    P2 = interp1d(setkin, np.dot(cvals,data2))(setkout) + b9/nd/km**2 * setkout**2 + b10/nd/km**2 * setkout**2
    P4 = interp1d(setkin, np.dot(cvals,data4))(setkout)
    
    return np.array([P0,P2,P4])
'''
  
def gelman_rubin_convergence(withinchainvar, meanchain, n, Nchains, ndim):
    
    """ Calculate Gelman & Rubin diagnostic
     1. Remove the first half of the current chains
     2. Calculate the within chain and between chain variances
     3. estimate your variance from the within chain and between chain variance
     4. Calculate the potential scale reduction parameter
   
    Inputs
    ------
        withinchainvar : array of the variances within each chains
        meanchain : array of the means within each chains
        n : length of the chains
        Nchains : number oc chains
        ndim : number of varied parameters
    Outputs
    ------
        The gelman rubin criteria
    
    
    """

    meanall  =  np.mean(meanchain, axis = 0)
    W  =  np.mean(withinchainvar, axis = 0)
    B  =  np.arange(ndim,dtype = np.float)
    for jj in range(0, ndim):
        B[jj]  =  0.
    for jj in range(0, Nchains):
        B  =  B + n*(meanall - meanchain[jj])**2/(Nchains-1.)
    estvar  =  (1. - 1./n)*W + B/n
    scalereduction  =  np.sqrt(estvar/W)

    return scalereduction
    
    
def match_para(theta, free_para, fix_para):
    
    """ Select the parameters that should be varied
    
    Inputs
    ------
    theta : array of all parameters
    free_para : list of boolean corresponding to the parameters to be varied
    fix_para : list of the values for the fixed parameters
    
    Outputs
    ------
    array of the parameters with the one not varied fixed to their fix_para value.
    
    """
    
    value_array  =  np.arange(len(free_para),dtype = np.float)
    
    counter  =  0
    for i in range(len(free_para)):
        if(free_para[i]  ==  True):
            value_array[i]  =  theta[counter]
            counter +=  1
        else: value_array[i]  =  fix_para[i]

    return value_array


def lnprior(theta, free_para, fix_para,bounds):
    
    """ Computes the prior 
    Inputs
    ------
    theta : array of all parameters
    free_para : list of boolean corresponding to the parameters to be varied
    fix_para : list of the values for the fixed parameters
    bounds : 2d array with the [min,max] values for the parameters
    
    Outputs
    ------
    the value of the log of the prior
    
    """
    
    value_array  =  match_para(theta, free_para, fix_para)
    Om = value_array[1]
    h = value_array[2]    
    withinprior = True
    for i in range(len(value_array)):
        withinprior = (withinprior) and (bounds[i][0] <= value_array[i] <= bounds[i][1])
        
    if withinprior:
        if withPlanck:
#            print(Om, h,  -0.5* rs(Om,h,f_fid))
            return -0.5* (RD-rs(Om,h,f_fid))**2/SIGMA_RD**2 
            
        else:            
            return 0.
    else:
     return -np.inf


def lnlike(theta,  kpred,chi2data,Cinvwdata,Cinvww, free_para, fix_para,bounds,Om_fid, sigma_prior = 100, marg_gaussian=False,binning=False,TableNkmu=None):
    
    """ Computes the log of the likelihood

    Inputs
    ------
    theta : array of all parameters
    kpred: the (logarithmically-spread) array of k on which Cinvww is evaluated
    chi2data : precomputed np.dot(Pdata,np.dot(Cinv,Pdata))
    Cinvwdata : precomputed np.dot(Cinvw,Pdata), where Cinvw is the inverse covariance multiplied by one window function
    Cinvww : the inverse of the covariance multiplied by the convolution matrices
    free_para : list of boolean corresponding to the parameters to be varied
    fix_para : list of the values for the fixed parameters
    bounds : 2d array with the [min,max] values for the parameters
    Om_fid : the value of cosmological parameter Om on the fiducial (for AP effect)
    interpolation_grid : the interpolation of the power spectra on the grid. Include the bispectrum if considered
    binning : whether to take into account the discreteness of the data power spectra. Must provide a number of modes per bin in TableNkmu
    TableNkmu : the number of modes per (k,mu) bin. Expected to be of the type kmean, mucentral, Nkmu
    Outputs
    ------
    the value of the log of the likelihood
    
    """
    
    # Because we have a precomputed grid for the cosmological parameters, need to check that we are within that grid (defined in the prior).
    if not np.isfinite(lnprior(theta, free_para, fix_para,bounds)):
        return -100000
    else :
    
        
        t0 = time.time()
        
        # If marginalization, no variation over b3,b5,b6,b7,b8,b9,b10,b11. 
        if marg_gaussian:
            #print(free_para)
            free_para = np.array(free_para)            
            free_para[[5,7,8,9,10,11,12,13]] = np.array([False]*8)
            fix_para[[5,7,8,9,10,11,12,13]] = np.array([0]*8)
            if withBisp:
                #b8 is not marginalized over when bisp is present
                free_para[[5,7,8,9,11,12,13]] = np.array([False]*7)
                fix_para[[5,7,8,9,11,12,13]] = np.array([0]*7)
        #print(match_para(theta, free_para, fix_para))
        lnAs,Om,h,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11 = match_para(theta, free_para, fix_para)
            
    # Import the power spectra interpolators on the grid
 
        if withBisp:    
            #Plininterp,Ploopinterp,Sigsqinterp,Bispinterp= interpolation_grid 
            Plininterp,Ploopinterp,Bispinterp= interpolation_grid 
        else: 
            #Plininterp,Ploopinterp,Sigsqinterp = interpolation_grid 
            Plininterp,Ploopinterp = interpolation_grid 

        
        kfull = Ploopinterp((lnAs,Om,h))[:,0]
    
        if check_if_multipoles_k_array(kfull):
            kfull = kfull[:len(kfull)/3] 
        kfullred = kfull[kfull<kpred.max()+0.1]
        Ploop = np.swapaxes(Ploopinterp((lnAs,Om,h)).reshape(3,len(kfull),22),axis1 = 1,axis2 = 2)[:,1:,:]
        Plin = np.swapaxes(Plininterp((lnAs,Om,h)).reshape(3,len(kfull),4),axis1 = 1,axis2 = 2)[:,1:,:]
        #sigsq = float(Sigsqinterp((lnAs,Om,h)))
        # Get sigma^2 to be removed
        
        if withBisp:
            if type(masktriangle) ==  type(None) or type(Bispdata) == type(None) or type(Bispinterp) ==  type(None):
                raise Exception('You want to use the bispectrum but forgot to provide a mask for the triangle or the data or the interpolation of the Bisp. Can be found in input/')
            if Cinv.shape[0] != xdata.shape + sum(masktriangle):
                raise Exception('You want to use the bispectrum but forgot to use the full covariance for power spectrum + Bisp. Can be found in input/Covariance')
        
            TermsBisp = Bispinterp((lnAs,Om,h))
            bval = np.array([1.,b1,b2,b4,b1*b11,b1**2,b1*b2,b1*b4,b1**3,b1**2*b2,b1**2*b4,b8**2])
            Bisp = 1.*np.dot(bval,TermsBisp[3:])
     
        # Compute the PS       
        valueb = np.array([b1,b2,b3,b4,b5,b6,b7,b8,b9,b10])        
        #Pmodel_original = computePS(valueb,Plin,Ploop,kfull,kfull,sigsq=sigsq)
        Pmodel_original = computePS(valueb,Plin,Ploop,kfull,kfull)
        Pmodel = Pmodel_original.copy()
    
        #The AP parameters
        qperp,qpar = get_AP_param(Om,Om_fid)
        
        if marg_gaussian:
            if withBisp:
                Pi_or = get_Pi_for_marg(Ploop,b1, TermsBisp[3:])
            else: 
                Pi_or = get_Pi_for_marg(Ploop,b1)
        #print(Pi_or.shape, 'is original Pi shape')
    
        if not binning:
            PmodelAP = APpowerspectraNkmu.changetoAPnobinning(Pmodel,kfull,kfullred,qperp,qpar,nbinsmu=100)
            if marg_gaussian:
                Pi_AP = APpowerspectraNkmu.changetoAPnobinningPi(Pi_or,kfull,kfullred,qperp,qpar,nbinsmu=100)

        else:
            if type(TableNkmu) == type(None):
                raise Exception('You want to account for binning but forgot to provide a TableNkmu (array of shape (3,n)) obtained from the sims/ Can be found in input/TableNkmu')
            else : 
                PmodelAP = APpowerspectraNkmu.changetoAPbinning(Pmodel,kfull,kfull,qperp,qpar,TableNkmu)
                if marg_gaussian:
                    Pi_AP = APpowerspectraNkmu.changetoAPbinningPi(Pi_or,kfull,kfull,qperp,qpar,TableNkmu)
        
        #print(kfullred.shape, kfull.shape, PmodelAP.shape)
        Pmodel_extrap = scipy.interpolate.interp1d(kfullred,PmodelAP,axis=-1,bounds_error=False,fill_value='extrapolate')(kpred)
        modelX = Pmodel_extrap.reshape(-1)
        if withBisp:
            modelX = np.concatenate([modelX, Bisp[masktriangle]])
        if marg_gaussian:
            Pi_extrap = (scipy.interpolate.interp1d(kfullred,Pi_AP,axis=-1,bounds_error=False,fill_value='extrapolate')(kpred)).reshape((Pi_AP.shape[0],-1))
            Pi_tot = 1.*Pi_extrap
            #print(Pi_tot.shape, ' is Pi total shape', kfull.shape, kfullred.shape, Pi_extrap.shape, ' is Pi extrap shape', kpred.shape)
            if withBisp:
                #Building a bigger Pi to include bispectrum term
                nparams = Pi_tot.shape[0]
                nkpred = Pi_tot.shape[1]
                nkbisp = sum(masktriangle)
                #print(Pi_tot.shape, TermsBisp[3:].shape, TermsBisp.shape)
                #Removing triangle contributions from the bispectrum expressions
                TermsBisp = TermsBisp[3:]
                #Applying mask and only getting b11 contribution
                bisp = TermsBisp[4][masktriangle]
                
                newPi = np.zeros(shape=(nparams+1, nkpred+nkbisp))
                newPi[:nparams, :nkpred] = Pi_tot
                newPi[-1, nkpred:] = b1*bisp
                #total Pi is now the correctly embedded one
                Pi_tot = 1.*newPi 

            Covbi = get_Covbi_for_marg(Pi_tot,Cinvww,sigma=sigma_prior)
            Cinvbi = np.linalg.inv(Covbi)
            vectorbi = np.dot(modelX,np.dot(Cinvww,Pi_tot.T))-np.dot(Cinvwdata,Pi_tot.T)
            chi2nomar = np.dot(modelX,np.dot(Cinvww,modelX))-2*np.dot(Cinvwdata,modelX)+chi2data
            chi2mar = -np.dot(vectorbi,np.dot(Cinvbi,vectorbi))+np.log(np.linalg.det(Covbi))
            chi2 = chi2mar + chi2nomar
        
        else:
             chi2 = np.dot(modelX,np.dot(Cinvww,modelX))-2*np.dot(Cinvwdata,modelX)+chi2data       
 
        if withPlanck:
            #print("Planck contrib is ", (rd-rs(Om,h,f_fid))**2/sigma_rd**2)
            return -0.5* ( chi2 )
        else: 
            return -0.5*chi2




def lnprob(theta, kpred,chi2data,Cinvwdata,Cinvww, free_para, fix_para,bounds,Om_fid, sigma_prior = 100, marg_gaussian=False,binning=False,TableNkmu=None):
   
    """ Computes the log of the probability (logprior + loglike)

    Inputs
    ------
    theta : array of all parameters
    kpred: the (logarithmically-spread) array of k on which Cinvww is evaluated
    chi2data : precomputed np.dot(Pdata,np.dot(Cinv,Pdata))
    Cinvwdata : precomputed np.dot(Cinvw,Pdata), where Cinvw is the inverse covariance multiplied by one window function
    Cinvww : the inverse of the covariance multiplied by the convolution matrices
    free_para : list of boolean corresponding to the parameters to be varied
    fix_para : list of the values for the fixed parameters
    bounds : 2d array with the [min,max] values for the parameters
    Om_fid : the value of cosmological parameter Om on the fiducial (for AP effect)
    interpolation_grid : the interpolation of the power spectra on the grid. Include the bispectrum if considered
    binning : whether to take into account the discreteness of the data power spectra. Must provide a number of modes per bin in TableNkmu
    TableNkmu : the number of modes per (k,mu) bin. Expected to be of the type kmean, mucentral, Nkmu
    Outputs
    ------
    the value of (logprior + loglike)
    
    """    
    lp  =  lnprior(theta, free_para, fix_para,bounds)
    
    if np.isfinite(lp) == False :
        dummy  =  -np.inf
        
    dummy  =  lp + lnlike(theta, kpred,chi2data,Cinvwdata,Cinvww, free_para, fix_para,bounds,Om_fid, sigma_prior=sigma_prior, marg_gaussian=marg_gaussian,binning=binning,TableNkmu=TableNkmu)

    return dummy


def lnlike2(theta,kpred,chi2data,Cinvwdata,Cinvww, free_para, fix_para, bounds,Om_fid, sigma_prior = 100, marg_gaussian=False,binning=False,TableNkmu=None):
    
    # Because we have a precomputed grid for the cosmological parameters, need to check that we are within that grid (defined in the prior).
    if not np.isfinite(1):
        return -100000
    else :
        #b3,b5,b6,b7,b8, b9,b10 = theta
        b3,b5,b6,b7,b8 = theta
        #b3,b5,b6,b7, b9,b10 = theta
        #b8 = 0
        lnAs,Om,h,b1,b2,b4,b9,b10= fix_para
            
    # Import the power spectra interpolators on the grid
 
        if withBisp:    
            #Plininterp,Ploopinterp,Sigsqinterp,Bispinterp= interpolation_grid 
            Plininterp,Ploopinterp,Bispinterp= interpolation_grid 
        else: 
            #Plininterp,Ploopinterp,Sigsqinterp = interpolation_grid 
            Plininterp,Ploopinterp = interpolation_grid 

        
        kfull = Ploopinterp((lnAs,Om,h))[:,0]
    
        if check_if_multipoles_k_array(kfull):
            kfull = kfull[:len(kfull)/3] 
        kfullred = kfull[kfull<kpred.max()+0.1]
        Ploop = np.swapaxes(Ploopinterp((lnAs,Om,h)).reshape(3,len(kfull),22),axis1 = 1,axis2 = 2)[:,1:,:]
        Plin = np.swapaxes(Plininterp((lnAs,Om,h)).reshape(3,len(kfull),4),axis1 = 1,axis2 = 2)[:,1:,:]
        #sigsq = float(Sigsqinterp((lnAs,Om,h)))
        # Get sigma^2 to be removed
        
        if withBisp:
            if type(masktriangle) ==  type(None) or type(Bispdata) == type(None) or type(Bispinterp) ==  type(None):
                raise Exception('You want to use the bispectrum but forgot to provide a mask for the triangle or the data or the interpolation of the Bisp. Can be found in input/')
            if Cinv.shape[0] != xdata.shape + sum(masktriangle):
                raise Exception('You want to use the bispectrum but forgot to use the full covariance for power spectrum + Bisp. Can be found in input/Covariance')
        
            TermsBisp = Bispinterp((lnAs,Om,h))
            bval = np.array([1.,b1,b2,b4,b1*b11,b1**2,b1*b2,b1*b4,b1**3,b1**2*b2,b1**2*b4,b8**2])
            Bisp = 1.*np.dot(bval,TermsBisp[3:])
     
        # Compute the PS       
        valueb = np.array([b1,b2,b3,b4,b5,b6,b7,b8,b9,b10])
        #Pmodel_original = computePS(valueb,Plin,Ploop,kfull,kfull,sigsq=sigsq)
        Pmodel_original = computePS(valueb,Plin,Ploop,kfull,kfull)
        Pmodel = Pmodel_original.copy()
    
        #The AP parameters
        qperp,qpar = get_AP_param(Om,Om_fid)
        
        if marg_gaussian:
            if withBisp:
                Pi_or = get_Pi_for_marg(Ploop,b1, TermsBisp[3:])
            else: 
                Pi_or = get_Pi_for_marg(Ploop,b1)
        #print(Pi_or.shape, 'is original Pi shape')
    
        if not binning:
            PmodelAP = APpowerspectraNkmu.changetoAPnobinning(Pmodel,kfull,kfullred,qperp,qpar,nbinsmu=100)
            if marg_gaussian:
                Pi_AP = APpowerspectraNkmu.changetoAPnobinningPi(Pi_or,kfull,kfullred,qperp,qpar,nbinsmu=100)

        else:
            if type(TableNkmu) == type(None):
                raise Exception('You want to account for binning but forgot to provide a TableNkmu (array of shape (3,n)) obtained from the sims/ Can be found in input/TableNkmu')
            else : 
                PmodelAP = APpowerspectraNkmu.changetoAPbinning(Pmodel,kfull,kfull,qperp,qpar,TableNkmu)
                if marg_gaussian:
                    Pi_AP = APpowerspectraNkmu.changetoAPbinningPi(Pi_or,kfull,kfull,qperp,qpar,TableNkmu)
        
        #print(kfullred.shape, kfull.shape, PmodelAP.shape)
        Pmodel_extrap = scipy.interpolate.interp1d(kfullred,PmodelAP,axis=-1,bounds_error=False,fill_value='extrapolate')(kpred)
        modelX = Pmodel_extrap.reshape(-1)
        if withBisp:
            modelX = np.concatenate([modelX, Bisp[masktriangle]])
        if marg_gaussian:
            Pi_extrap = (scipy.interpolate.interp1d(kfullred,Pi_AP,axis=-1,bounds_error=False,fill_value='extrapolate')(kpred)).reshape((Pi_AP.shape[0],-1))
            Pi_tot = 1.*Pi_extrap
            #print(Pi_tot.shape, ' is Pi total shape', kfull.shape, kfullred.shape, Pi_extrap.shape, ' is Pi extrap shape', kpred.shape)
            if withBisp:
                #Building a bigger Pi to include bispectrum term
                nparams = Pi_tot.shape[0]
                nkpred = Pi_tot.shape[1]
                nkbisp = sum(masktriangle)
                #print(Pi_tot.shape, TermsBisp[3:].shape, TermsBisp.shape)
                #Removing triangle contributions from the bispectrum expressions
                TermsBisp = TermsBisp[3:]
                #Applying mask and only getting b11 contribution
                bisp = TermsBisp[4][masktriangle]
                
                newPi = np.zeros(shape=(nparams+1, nkpred+nkbisp))
                newPi[:nparams, :nkpred] = Pi_tot
                newPi[-1, nkpred:] = b1*bisp
                #total Pi is now the correctly embedded one
                Pi_tot = 1.*newPi 

            Covbi = get_Covbi_for_marg(Pi_tot,Cinvww,sigma=sigma_prior)
            Cinvbi = np.linalg.inv(Covbi)
            vectorbi = np.dot(modelX,np.dot(Cinvww,Pi_tot.T))-np.dot(Cinvwdata,Pi_tot.T)
            chi2nomar = np.dot(modelX,np.dot(Cinvww,modelX))-2*np.dot(Cinvwdata,modelX)+chi2data
            chi2mar = -np.dot(vectorbi,np.dot(Cinvbi,vectorbi))#+np.log(np.linalg.det(Covbi))
            chi2 = chi2mar + chi2nomar
        
        else:
             chi2 = np.dot(modelX,np.dot(Cinvww,modelX))-2*np.dot(Cinvwdata,modelX)+chi2data       
 
        if withPlanck:
            #print("Planck contrib is ", (rd-rs(Om,h,f_fid))**2/sigma_rd**2)
            return -0.5* ( chi2 )
        else: 
            return -0.5*chi2
                                                                                                    ###########################################
                                                                                                    ###  Main program  ########################
                                                                                                    ###########################################


if __name__ ==  "__main__":

    # Table of cosmological parameters according to seems
    print("started")
    dfcosmo = pd.read_csv(opa.join(INPATH,'DataFrameCosmosims.csv'),index_col=0)
    boxnumber = sys.argv[1]
    KMAX = float(sys.argv[2])
    #simtype = "LightConeDida"
    simtype = sys.argv[3]

    planckchain = int(sys.argv[4])
    
    series_cosmo = dfcosmo.loc[simtype]
    if planckchain == 1:
        withPlanck = True    
    else:
        withPlanck = False
    #try:
    #    sys.argv[5]
    # except:
    #    raise Exception("You asked for a non-default grid but didn't give its name as an argument!")
    gridname = str(sys.argv[5]) 
    print(gridname)
    # Load the row that we are interested in
    
    
    # gridname = series_cosmo.loc['gridname']#
    
    

    # COSMOLOGICAL GLOBALS: fiducial model (should match input sim data!)
    Om_fid  =  series_cosmo.loc['Omega_m']
    lnAs_fid = series_cosmo.loc['lnAs']
    h_fid  =  series_cosmo.loc['h']
    #z_pk = 0.57
    z_pk = series_cosmo.loc['z_pk']

    

    ob_fid = series_cosmo.loc['omega_b']
    # ratio omega_b/omega_c
    f_fid = ob_fid / (Om_fid*h_fid**2 - ob_fid)

    #withPlanck = False 


    
    
    withBisp = False 
    
    window = True
    #### Choice for the data #####
    #For lightcone simulations, need to specify north or south for now (later, merge the two but I'm missing the covariance for SGC
    #Change this back when not doing Challenge boxes
    ZONE = 'NGC'
    if "Challenge" in simtype:
        ZONE = ''
    if "Hector" not in simtype:
        withBisp = False

    kmin = 0.01
    kminbisp = kmin
    
    kmaxbisp = float(sys.argv[6])
    if kmaxbisp > 0:
        withBisp = True
        print('kmax bispectrum is bigger than zero, withBisp is %s'%withBisp)
    #withMarg is 1 or 0
    withMarg = int(sys.argv[7])

    #Defines if will use a covariance with theory error or not
    theoryCovness = float(sys.argv[8])

    #Priors max
    priorsup = float(sys.argv[9])

    simtype2 = simtype+ZONE
    runtype = simtype+ZONE+'prior%s'%priorsup
    #workaround setting of marg_gauss, kmaxbisp = 0.07 is true and kmaxbisp = 0.08 is false
    if withMarg:
        marg_gaussian = True
        #b8 is free if bispectrum, b11 is never varied now (marged if bisp)
        free_para =  [True,True,True,True,True,False,True,False,False,False,withBisp,False,False,False]
        a = 1.8
        print("withMarg is ", withMarg, " setting marg_gaussian to ", marg_gaussian)
    else:
        marg_gaussian = False
        free_para =  [True,True,True,True,True,True,True,True,True,True,False,True,False,withBisp]
        a = 1.15
        print("withMarg is", withMarg, " setting marg_gaussian to ", marg_gaussian)
    if ZONE != '':    
        dataQ = np.loadtxt(opa.join(INPATH,'Window_functions/dataQ_%s.txt'%ZONE)).T 
    elif 'ChallengeQuarter' in simtype:
        dataQ = np.loadtxt(opa.join(INPATH,'Window_functions/dataQ_ChallengeQuarter.dat')).T
    if 'Challenge' in simtype:
        #CHANGE THIS BACK, SETTING FULL COV TO QUARTER
        #change RD to match the CHallenge Box cosmology

        RD = 147.68188
        if 'Quarter' and 'Japan' not in simtype:
            window = False
            simtype_false = 'ChallengeQuarter'+boxnumber
            print('Using quarter covariance instead of full')
            Full_Cov = np.loadtxt(opa.join(INPATH,'Covariance/Cov%s%s.dat'%(simtype_false,ZONE)))
        elif 'Japan' in simtype:
            Full_Cov = np.loadtxt(opa.join(INPATH,'Covariance/Cov%s_%s.dat'%(simtype,boxnumber)))
        else:
            Full_Cov = np.loadtxt(opa.join(INPATH,'Covariance/Cov%s%s.dat'%(simtype,ZONE)))
    else:
        if theoryCovness == 0:
            Full_Cov = np.loadtxt(opa.join(INPATH,'Covariance/Cov%s%sdata.dat'%(simtype,ZONE)))
        else:
            Full_Cov = np.loadtxt(opa.join(INPATH,'Covariance/Cov%s%s_Theory_a%s.dat'%(simtype,ZONE,theoryCovness)))
            runtype+='tCov%s'%theoryCovness
    print("here")    
    #bra
    if not marg_gaussian:
        if withBisp:
            raise(Exception("Non-marginalized bispectrum is currently not implemented!"))
    Bispdata = [] 
    masktriangle = []
    if withBisp:
        runtype += 'withBispkmax%s'%kmaxbisp
        Full_Cov = np.loadtxt(opa.join(INPATH,'Covariance/Cov%s%s_Bisp.dat'%(simtype,ZONE)))   
        Q1,Q2,Q3,Bispdata = np.loadtxt(opa.join(INPATH,'DataSims/Bispred_LightConeHector_%s_%s.dat'%(ZONE,boxnumber))).T
        masktriangle = (Q1>=kminbisp)&(Q1<=kmaxbisp)&(Q1<=Q2+Q3)&(Q1>=abs(Q2-Q3))&(Q2>=kminbisp)&(Q2<=kmaxbisp)&(Q3>=kminbisp)&(Q3<=kmaxbisp)
    #print(masktriangle.shape)

                                                    
                                                    ####### OPTIONS FOR MCMC #########
    
    binning = False
    TableNkmu = None
    #marg_gaussian = True
    
    nbinsAs = 100
    nbinsOm = 48
    nbinsh = 48

    if 'Patchy' in gridname:
        nbinsAs = 70
        nbinsOm = 48
        nbinsh = 72
    

    lnAsmin,lnAsmax,Ommin,Ommax,hmin,hmax,interpolation_grid = get_grid(gridname,nbinsAs=nbinsAs,nbinsOm=nbinsOm,nbinsh=nbinsh,withBisp=withBisp)  
    print("got grid!")     
##############################
###  Priors ###################
#############################

    #### The uniform prior on the b_i#####
    bmin = -priorsup
    bmax = priorsup

    # We require b_1>0 and b_4 is large due to the PCA that we do between b2 and b4
    bmintab = [0, bmin, bmin, -100, bmin, bmin, bmin, -100, bmin, bmin, -100]
    bmaxtab = [bmax, bmax, bmax, 100, bmax, bmax, bmax, 100, bmax, bmax, 100]
    

    ##### The bounds for the minimizer ######
    bfmaxtab = np.concatenate([[lnAsmax,Ommax,hmax],bmaxtab])
    bfmintab = np.concatenate([[lnAsmin,Ommin,hmin],bmintab])

    bounds = zip(bfmintab,bfmaxtab)

    ##### Initial guess for the b_i #####
    '''
    inipos = np.array([1.85 ,  -2.62623719,  -0.39661384,   4.21514113, 
    8.36786486, -29.68630616,   1.03528956, 
    0, 40.00717862,  0,   100])
    '''
    inipos = np.array([ 1.85 ,  1,  1,  1, 
                        1, 1,   1, 
                        0, 1,  1,   10])

    ##### Guess for the \sigma, to help with the initial position of walkers #####
    onesigma = np.array([   2.48736140e-01,   4.40317511e-02,   2.65820186e-02, 
                            0.1,   0.1,   0.1, 0.1,  
                            0.1,   0.1, 0.1,  
                            1,   0.1, 0.1, 1]) 

    
    kmaxtab = [KMAX]

    print("Starting process")
    print("bounds are ")
    print (bounds)
    sys.stdout.flush()
    kmaxname = ['kmax%s'%kmax for kmax in kmaxtab]

    for boxnumber in [boxnumber]:
        if 'Challenge' in simtype:
            if 'Japan' not in simtype:
                simtype = simtype[:-1]
        kPS,PSdata,_ = np.loadtxt(opa.join(INPATH,'DataSims/ps1D_%s%s_%s.dat'%(simtype,ZONE,boxnumber))).T
        klog = np.loadtxt(opa.join(INPATH,'Window_functions/k.dat'))
        for indexk,kmax in enumerate(kmaxtab):    

            xdata = kPS[(kPS<kmax)&(kPS>kmin)]
            ydata = PSdata[(kPS<kmax)&(kPS>kmin)]
            indexkred = np.argwhere((kPS<kmax)&(kPS>kmin))[:,0]
            
            kindex = 1.*indexkred
            kmask = np.array([False]*len(kPS))
            kmask[kindex.astype(int)] = True
            kmask = kmask[:len(kmask)/3]

            if withBisp:
                indextriangle = np.argwhere(masktriangle)[:,0]+kPS.shape[0]
                indexkred = np.concatenate([indexkred,indextriangle])
                ydata = np.concatenate([ydata, Bispdata[masktriangle]])
            print("triangle length is", np.argwhere(masktriangle)[:,0].shape)        
            Covred = Full_Cov[indexkred[:,None],indexkred]

            Cinv = np.linalg.inv(Covred)
            #Window function shouldn't be applied for full challenge box
            if 'Challenge' in simtype and 'Quarter' not in simtype:
                #print("Am I doing this instead?")
                kpred = xdata[:len(xdata)/3]
                Cinvw = Cinv
                Cinvww = Cinv
            elif 'Quarter' in simtype:
                raise(Exception("Quarter window function not yet implemented! Look at hard-coded path in WindowFunctionFourier.py for more information"))
            else:
                #print("I'm applying the window function!!")
                kpred,Cinvw,Cinvww = WindowFunctionFourier.apply_window_covariance(Cinv,xdata,thin=2, bisp = withBisp, indexkred = kmask, masktriangle = masktriangle)
            
            chi2data = np.dot(ydata,np.dot(Cinv,ydata))
            Cinvwdata = np.dot(ydata, Cinvw)     

        

    ##################################################################
    ## Setting up second fit ###########
    ##################################################################
    
    lnlike = np.load('/exports/pierre/EFTofBOSS/output/lnlikechain'+simtype2+'prior'+str(priorsup)+'gaussMargbox_'+str(boxnumber)+'kmax_'+str(kmax)+'run_0.npy')
    chain = np.array(np.load('/exports/pierre/EFTofBOSS/output/samplerchain'+simtype2+'prior'+str(priorsup)+'gaussMargbox_'+str(boxnumber)+'kmax_'+str(kmax)+'run_0.npy'))
    id1,idx=unravel_index(lnlike.argmax(), lnlike.shape)
    print (lnlike.argmax())
    id2,idy=unravel_index(lnlike.argmin(), lnlike.shape)
    print (chain[id2,idy*5,:])
    print (chain[id1,idx*5,:])
    print ('fid As%s, Om%s, h%s'%(lnAs_fid,Om_fid,h_fid))
    lnAs,Om,h,b1,b2,b4 = chain[id1,idx*5,:]
    b9 = 0
    b10 = 0

    free_ml = np.array([lnAs,Om,h,b1,b2,b4,0,0])

    free_para2 =  [False,False,False,False,False,True,False,True,True,True,True,False,False,False]
    
    # b8 = 0
    #free_para2 =  [False,False,False,False,False,True,False,True,True,True,False,True,True,False]

    # non-margi
    #free_para2 = [True,True,True,True,True,True,True,True,True,True,False,True,True,False]
    nparam2 = len(free_para2)
    

    all_true2  =  np.concatenate(([lnAs,Om,h,b1,b2,1, b4], 7*[0] ))
    #all_name2  =  np.concatenate(([r'$A_s$',r'$\Omega_m$',r'$h$'],[r'$b_%s$'%(i+1) for i in range(len(inipos))]))
    
    ndim2  =  sum(free_para2)
    #print(ndim)
    #print(free_para2)
    #print(fix_para2)
    #free_true2 = all_true2[free_para2]
                #################################
    ## Find maximum likelihood ######
    #################################

    chi2bis  =  lambda theta2: -2 * lnlike2(theta2, kpred,chi2data,Cinvwdata,Cinvww, free_para2, free_ml, bounds, Om_fid, sigma_prior = priorsup, binning=binning,marg_gaussian=False,TableNkmu=TableNkmu)

    result2  =  op.minimize(chi2bis, np.array(all_true2)[free_para2], method = 'SLSQP',bounds = np.array(bounds)[free_para2],options = {'maxiter':100})
    
    all_ml2  =  result2["x"]
    
    free_ml2 = all_ml2[:ndim2]

    minchi22  =  result2["fun"]
         
    print(result2)
    #print("ap params", get_AP_param(free_ml[1], Om_fid))
    if type(masktriangle) == type(None):
        dof2 = len(xdata) - ndim2
    else:
        dof2 = len(xdata) + sum(masktriangle) - ndim2
    print('minchi2 = ' + str(minchi22), dof2)
    #np.savetxt(opa.join(OUTPATH,"minchi2%sbox_%skmax_%s.txt")%(runtype,boxnumber,kmax),np.concatenate([free_ml,[minchi2,dof]]))

    #################################
    ## Produce spectra ######
    #################################

    b3,b5,b6,b7,b8 = free_ml2
    #b3,b5,b6,b7, b8, b9,b10 = free_ml2
    #b8 = 0

    #lnAs,Om,h,b1,b2,b3,b4,b5,b6,b7,b9,b10 = free_ml2


    Plininterp,Ploopinterp = interpolation_grid 
        
    kfull = Ploopinterp((lnAs,Om,h))[:,0]
    
    if check_if_multipoles_k_array(kfull):
        kfull = kfull[:len(kfull)/3] 
    kfullred = kfull[kfull<kpred.max()+0.1]
    Ploop = np.swapaxes(Ploopinterp((lnAs,Om,h)).reshape(3,len(kfull),22),axis1 = 1,axis2 = 2)[:,1:,:]
    Plin = np.swapaxes(Plininterp((lnAs,Om,h)).reshape(3,len(kfull),4),axis1 = 1,axis2 = 2)[:,1:,:]
        
    valueb = np.array([b1,b2,b3,b4,b5,b6,b7,b8,b9,b10])        
    Pmodel = computePS(valueb,Plin,Ploop,kfull,kfull)

    qperp,qpar = get_AP_param(Om,Om_fid)

    #qperp = 1
    #qpar =1

    PmodelAP = APpowerspectraNkmu.changetoAPnobinning(Pmodel,kfull,kfullred,qperp,qpar,nbinsmu=100)
    if 'LightConeHector' in simtype2:
        Pmodel_extrap = WindowFunctionFourier.apply_window_PS(kfullred,PmodelAP,xdata[:len(xdata)/3],withmask=True,windowk=0.1)
        modelX = Pmodel_extrap.reshape(-1)
    else:
        Pmodel_extrap = scipy.interpolate.interp1d(kfullred,PmodelAP,axis=-1,bounds_error=False,fill_value='extrapolate')(xdata[:len(xdata)/3])
        modelX = Pmodel_extrap.reshape(-1)

    np.savetxt(opa.join(OUTPATH,"fit%sbox_%skmax_%s.txt")%(runtype,boxnumber,kmax), zip(xdata, ydata, modelX, np.diag(Covred)))
    np.savetxt(opa.join(OUTPATH,"bestfit%sbox_%skmax_%s.txt")%(runtype,boxnumber,kmax), [lnAs,Om,h,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10])
    np.savetxt(opa.join(OUTPATH,"minchi2%sbox_%skmax_%s.txt")%(runtype,boxnumber,kmax),[minchi22,dof2])