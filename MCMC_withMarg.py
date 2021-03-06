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
import WindowFFTlog

# Data paths
INPATH = opa.abspath(opa.join(THIS_PATH,'input'))
#INPATH2 = opa.abspath('/scratch/users/kokron/')
OUTPATH = opa.abspath(opa.join(THIS_PATH,'output')) 

###########################################
###  Functions  ###########################
###########################################


#####################################################################################################################################################
################## Sound horizon at decoupling #######################
###
### Option 'withPlanck' added (l.433) in lnlike: setup it up at l.518.
###
################################################### - Pierre, 10/08/18
#####################################################################################################################################################
# speed of light [km/s]
c = 299792.458 
# omega_gamma = Omega_gamma h^2: normalized physical photon density today (T_cmb = 2.7255 (CLASS))
og = 2.47282e-5
#Nur: Number of ultra-relativistic species (CLASS):
Nur = 3.046 
# omega_radiation
orad = (1.+ Nur*7./8.*(4./11.)**(4./3.))*og
# Baryon-photon decoupling redshift (PLANCK 2015 TT,TE,EE+lowP+lensing (Table 4)):
zd = 1059.62
sigma_zd = 0.31 # rd(zd+sigma)-dr(zd-sigma) < 0.2 sigma_rd: we take zd to be a delta function
# Sound horizon at decoupling [Mpc] (PLANCK 2015 TT,TE,EE+lowP+lensing (Table 4)):
rd = 147.41
sigma_rd = 0.30

def rs(Om,h,f_fid):
    om = Om*h**2
    ob = om * f_fid/(f_fid+1.) # f_fid: fiducial ratio omega_b/omega_c
    R = 0.75 * ob/og
    return 2.*c/100./np.sqrt(3.*R*om) * np.log( ( np.sqrt(1.+zd+R) + np.sqrt((1.+zd)*R*orad/om+R) ) / np.sqrt(1.+zd) / (1.+np.sqrt(R*orad/om)) ) ;
#####################################################################################################################################################


def get_Pi_for_marg(Ploop,b1):
    Pi = np.array([Ploop[:,3,:]+b1*Ploop[:,7,:],
             Ploop[:,15,:]+b1*Ploop[:,12,:],
             Ploop[:,16,:]+b1*Ploop[:,13,:],
             Ploop[:,17,:]+b1*Ploop[:,14,:],
             Ploop[:,18,:],
             Ploop[:,19,:],
             Ploop[:,20,:]])
    
    
    
    return Pi
    

def get_Covbi_for_marg(Pi_data,Cinv,sigma=200):

    Covbi = np.dot(Pi_data,np.dot(Cinv,Pi_data.T))+ 1./sigma**2*np.identity(Pi_data.shape[0])
    return Covbi
    
    

def apply_WF_bi(Pi_data_AP,setkin,setkout,dataQ,kr=0.5,damp=True,extrap=True,k_junc_low=0.01,k_junc_high=0.2,ktr=100,sig=0.5,a=2,dlnx=0.1):
    
    if check_if_multipoles_k_array(setkin):
        setkin = setkin[:len(setkin)/3]
    
    Pi_data_final_nob8 = np.array([WindowFFTlog.transformQ(np.concatenate(PAP),setkin,setkout,dataQ,kr=kr,damp=damp,
                                                        extrap=extrap,k_junc_low=k_junc_low,k_junc_high=k_junc_high,ktr=ktr,sig=sig,a=a,dlnx=dlnx,BBKSlow=False,nmax=10) for PAP in Pi_data_AP[:-3]])
    Pi_data_final_nob8 = Pi_data_final_nob8.reshape((4,-1))
    
    Pi_data_final_b8 = scipy.interpolate.interp1d(setkin,Pi_data_AP[-3:],axis=-1)(setkout[:len(setkout)/3]).reshape((3,-1))
    Pi_data_final = np.concatenate([Pi_data_final_nob8,Pi_data_final_b8])
    return Pi_data_final


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
            


def get_grid(gridname,nbinsAs=100,nbins = 50,withBisp=False):
    
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

    theta3D = thetatab.reshape((nbinsAs,nbins,nbins,3))

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
    Tablesigsq = np.load(opa.abspath(opa.join(INPATH,'GridsEFT/Tablesigsq%s.npy'%gridname)))

    Plininterp = scipy.interpolate.RegularGridInterpolator((lnAstab,Omtab,htab),TablePlin.reshape((nbinsAs,nbins,nbins,TablePlin.shape[-2],TablePlin.shape[-1])))
    Ploopinterp = scipy.interpolate.RegularGridInterpolator((lnAstab,Omtab,htab),TablePloop.reshape((nbinsAs,nbins,nbins,TablePloop.shape[-2],TablePloop.shape[-1])))
    Sigsqinterp = scipy.interpolate.RegularGridInterpolator((lnAstab,Omtab,htab),Tablesigsq.reshape((nbinsAs,nbins,nbins)))
    
    interpolations = [Plininterp,Ploopinterp,Sigsqinterp]
    if withBisp:
        TableBisp = np.load(opa.abspath(opa.join(INPATH,'GridsEFT/TableBisp%s.npy'%gridname)))
        Bispinterp = scipy.interpolate.RegularGridInterpolator((lnAstab,Omtab,htab),TableBisp.reshape((nbinsAs,nbins,nbins,TableBisp.shape[-2],TableBisp.shape[-1])))
        interpolations = [Plininterp,Ploopinterp,Sigsqinterp,Bispinterp]
        
    return lnAsmin,lnAsmax,Ommin,Ommax,hmin,hmax,interpolations
    
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
    b1,b2,b3,b4,b5,b6,b7,b8,b9,b10 = cvals#
    
    # the columns of the Ploop data files.
    cvals = np.array([1,b1,b2,b3,b4,b1*b1,b1*b2,b1*b3,b1*b4,b2*b2,b2*b4,b4*b4,b1*b5,b1*b6,b1*b7,b5,b6 ,b7,b8 ,b9,b10])
    
    # Check the k-arrays are in the right format (not concatenated for multipoles)
    if check_if_multipoles_k_array(setkin):
        setkin = setkin[:len(setkin)/3]
    if check_if_multipoles_k_array(setkout):
        setkout = setkin[:len(setkout)/3]
         
        
    P0 = interp1d(setkin,np.dot(cvals,data0)+datalin0[0]+b1*datalin0[1]+b1*b1*datalin0[2] - 2*(-b1 + b2 + b4)**2*sigsq)(setkout)
    P2 = interp1d(setkin,np.dot(cvals,data2)+datalin2[0]+b1*datalin2[1]+b1*b1*datalin2[2])(setkout)
    P4 = interp1d(setkin,np.dot(cvals,data4)+datalin4[0]+b1*datalin4[1]+b1*b1*datalin4[2])(setkout)
    
    return np.array([P0,P2,P4])


  
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
    
    withinprior = True
    for i in range(len(value_array)):
        withinprior = (withinprior) and (bounds[i][0] <= value_array[i] <= bounds[i][1])
        
    if withinprior:            
        return 0.
    else:
     return -np.inf


def lnlike(theta, xdata, ydata, Cinv, free_para, fix_para,bounds,Om_fid, f_fid,marg_gaussian=False,binning=False,TableNkmu=None, window=True,dataQ=None,withBisp=False,masktriangle=None,Bispdata=None,withPlanck=False):
    
    """ Computes the log of the likelihood
    Inputs
    ------
    theta : array of all parameters
    xdata : the k on which the data is evaluated (concatenated for multipoles)
    ydata : the data concatenated for multipoles
    Cinv : the inverse of the covariance
    free_para : list of boolean corresponding to the parameters to be varied
    fix_para : list of the values for the fixed parameters
    bounds : 2d array with the [min,max] values for the parameters
    fiducial : the value of cosmological parameters on the fiducial (for AP effect)
    interpolation_grid : the interpolation of the power spectra on the grid. Include the bispectrum if considered
    binning : whether to take into account the discreteness of the data power spectra. Must provide a number of modes per bin in TableNkmu
    TableNkmu : the number of modes per (k,mu) bin. Expected to be of the type kmean, mucentral, Nkmu
    window : whether to take into account the window function (when not a periodic square box, such as Lightcone). Must provide the data for the Q matrices (see e.g. section 4.3.2 of arXiv:1706.02362)
    dataQ : the Q matrices (multipoles of the window function) taken from Florian. Of the form stab,dummy,Q0,Q2,Q4,Q6,Q8,mueff
    withBisp : whether we also include the Bisepctrum. In this case one need to provide a mask for the triangle we consider, as well as the data.
    masktriangle : mask for the triangle to consider
    Bispdata : the data for the bispectrum
    Outputs
    ------
    the value of the log of the likelihood
    
    """
    
    # Because we have a precomputed grid for the cosmological parameters, need to check that we are within that grid (defined in the prior).
    if not np.isfinite(lnprior(theta, free_para, fix_para,bounds)):
        return -100000
    else :
        t0 =  time.time()
        lnAs,Om,h,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11 = match_para(theta, free_para, fix_para)
        
    # Import the power spectra interpolators on the grid
        if withBisp:
            Plininterp,Ploopinterp,Sigsqinterp,Bispinterp = interpolation_grid
        else:
            Plininterp,Ploopinterp,Sigsqinterp = interpolation_grid 
            Bispinterp = None
        
        kfull = Ploopinterp((lnAs,Om,h))[:,0]
    
        if check_if_multipoles_k_array(kfull):
            kfull = kfull[:len(kfull)/3] 
        
        Ploop = np.swapaxes(Ploopinterp((lnAs,Om,h)).reshape(3,len(kfull),22),axis1 = 1,axis2 = 2)[:,1:,:]
        Plin = np.swapaxes(Plininterp((lnAs,Om,h)).reshape(3,len(kfull),4),axis1 = 1,axis2 = 2)[:,1:,:]         
        
        # Get sigma^2 to be removed
        sigsq = Sigsqinterp((lnAs,Om,h))  
        
        
        
     
        # Compute the PS       
        valueb = np.array([b1,b2,b3,b4,b5,b6,b7,b8,b9,b10])        
        Pmodel_original = computePS(valueb,Plin,Ploop,kfull,kfull,sigsq=sigsq)
        Pmodel = Pmodel_original.copy()
    
        #The AP parameters
        qperp,qpar = get_AP_param(Om,Om_fid)
        
        if marg_gaussian:
            Pi_data_or = get_Pi_for_marg(Ploop,b1)
    
        if not binning:
            PmodelAP = APpowerspectraNkmu.changetoAPnobinning(Pmodel,kfull,kfull,qperp,qpar)
            Pi_data_AP = np.array([APpowerspectraNkmu.changetoAPnobinning(P,kfull,kfull,qperp,qpar) for P in Pi_data_or])
        else:
            if type(TableNkmu) == type(None):
                raise Exception('You want to account for binning but forgot to provide a TableNkmu (array of shape (3,n)) obtained from the sims/ Can be found in input/TableNkmu')
            else : 
                PmodelAP = APpowerspectraNkmu.changetoAPbinning(Pmodel,kfull,kfull,qperp,qpar,TableNkmu)
                Pi_data_AP = np.array([APpowerspectraNkmu.changetoAPbinning(P,kfull,kfull,qperp,qpar,TableNkmu) for P in Pi_data_or])
            
    
    
        if window:
            if type(dataQ) ==  type(None):
                raise Exception('You want to account for window function but forgot to provide a dataQ (array of shape (8,n)) obtained from the sims. Can be found in input/dataQ')
            else: 
                k_junc_high = 0.27
                k_junc_low = kfull[1]
                
                # In principle there is a damping function applied to the hole power spectrum. It is a step function controlled by
                # a k transition ktr and a slope sigma. With the choice below, the function as no effect.
                ktr = 100
                sig = 0.5
                
                #This is for the FFT (essentially, the product of the midpoints of the k and r arrays)
                kr = 0.5
                
                #This ensures that at high k, the extrapolation is of the for b(k/k_high)^c x exp(-a(k-khigh)/khigh), where a is positive and damps the power spectrum
                
                damp = True
                            
                Pmodelfinal = WindowFFTlog.transformQ(np.concatenate(PmodelAP),np.concatenate([kfull,kfull,kfull]),xdata,dataQ,kr=kr,damp=damp,
                                                        extrap=True,k_junc_low=k_junc_low,k_junc_high=k_junc_high,ktr=ktr,sig=sig,a=2)
                if marg_gaussian:
                    Pi_data_final = apply_WF_bi(Pi_data_AP,np.concatenate([kfull,kfull,kfull]),xdata,dataQ,kr=kr,damp=damp,
                                                        extrap=True,k_junc_low=1.1*k_junc_low,k_junc_high=k_junc_high,ktr=ktr,sig=sig,a=10,dlnx=0.1)
        else:
            Pmodelfinal = APpowerspectraNkmu.changetoAPnobinning(PmodelAP,kfull,xdata,1,1) #This is just to interpolate the power spectrum on xdata
            if marg_gaussian:
                Pi_data_final = np.array([APpowerspectraNkmu.changetoAPnobinning(PAP,kfull,xdata,1,1) for PAP in Pi_data_AP])
                Pi_data_final = Pi_data_final.reshape((Pi_data_final.shape[0],-1))
        
        
        modelX = np.concatenate(Pmodelfinal)
                    
        if withBisp:
            if type(masktriangle) ==  type(None) or type(Bispdata) == type(None) or type(Bispinterp) ==  type(None):
                raise Exception('You want to use the bispectrum but forgot to provide a mask for the triangle or the data or the interpolation of the Bisp. Can be found in input/')
            if Cinv.shape[0] != xdata.shape + sum(masktriangle):
                raise Exception('You want to use the bispectrum but forgot to use the full covariance for power spectrum + Bisp. Can be found in input/Covariance')
        
            TermsBisp = Bispinterp((lnAs,Om,h))
            bval = np.array([1.,b1,b2,b4,b1*b11,b1**2,b1*b2,b1*b4,b1**3,b1**2*b2,b1**2*b4,b8**2])
            Bisp = 1.*np.dot(bval,TermsBisp[3:])
        
            modelX = np.concatenate([modelX,Bisp[masktriangle]])
            ydata = np.concatenate([ydata,Bispdata[masktriangle]])
        
        
        
        
    
        diff  =  (modelX - ydata)
        step1 = np.dot(Cinv,diff)
        chi2 = np.dot(diff,step1)
        
        
        if marg_gaussian:
            Covbi = get_Covbi_for_marg(Pi_data_final,Cinv)
            Cinvbi = np.linalg.inv(Covbi)
            vectorbi = np.dot(Pi_data_final,step1)
            
            chi2mar = -np.dot(vectorbi,np.dot(Cinvbi,vectorbi))+np.log(np.linalg.det(Covbi))
            
            chi2 = chi2 + chi2mar
        
        if np.isnan(chi2) or chi2>1000:

            modelX_original = np.concatenate(scipy.interpolate.interp1d(kfull,PmodelAP,axis=-1)(xdata[:len(xdata)/3]))#np.concatenate(APpowerspectraNkmu.changetoAPnobinning(Pmodel_original,kfull,xdata,qperp,qpar))
 
            if withBisp:
                modelX_original = np.concatenate([modelX_original,Bisp[masktriangle]])
            diff_original  =  (modelX_original - ydata)
            step1_original = np.dot(Cinv,diff_original)
            chi2_original = np.dot(diff_original,step1_original)
            
            if marg_gaussian:
                Pi_data_original = scipy.interpolate.interp1d(kfull,Pi_data_AP,axis=-1)(xdata[:len(xdata)/3]).reshape((Pi_data_AP.shape[0],-1))
                Covbior = get_Covbi_for_marg(Pi_data_original,Cinv)
                Cinvbior = np.linalg.inv(Covbior)
                
                vectorbior = np.dot(Pi_data_original,step1_original)  
                    
                        
                chi2_original = np.dot(diff_original,step1_original) - np.dot(vectorbior,np.dot(Cinvbior,vectorbior))+np.log(np.linalg.det(Covbior))

            if chi2_original<200:
                
                ntry = 0
                nmax = 5
                khighlist = np.concatenate([np.linspace(0.9*k_junc_high,k_junc_high,nmax),np.linspace(k_junc_high,1.1*k_junc_high,nmax+1)[1:]])
                while ((np.isnan(chi2) or abs(chi2_original-chi2)>100) and ntry<10):
                    
                    k_junc_high =  khighlist[nmax-1+(-1)**(ntry+1)*(ntry+1)/2]

                
                    Pmodelfinal = WindowFFTlog.transformQ(np.concatenate(PmodelAP),np.concatenate([kfull,kfull,kfull]),xdata,dataQ,kr=kr,damp=damp,extrap=True,k_junc_low=k_junc_low,k_junc_high=k_junc_high,ktr=ktr,sig=sig,a=2)
                    modelX = np.concatenate(Pmodelfinal)
                    if withBisp:
                        modelX = np.concatenate([modelX,Bisp[masktriangle]])
                    diff  =  (modelX - ydata)
                    step1 = np.dot(Cinv,diff)
                    chi2 = np.dot(diff,step1)
                    ntry += 1
                    if marg_gaussian:
                        Pi_data_final = apply_WF_bi(Pi_data_AP,np.concatenate([kfull,kfull,kfull]),xdata,dataQ,kr=kr,damp=damp,
                                                        extrap=True,k_junc_low=1.1*k_junc_low,k_junc_high=k_junc_high,ktr=ktr,sig=sig,a=10,dlnx=0.1)
                        Covbi = get_Covbi_for_marg(Pi_data_final,Cinv)
                        Cinvbi = np.linalg.inv(Covbi)
                        vectorbi = np.dot(Pi_data_final,step1)
                        chi2mar = -np.dot(vectorbi,np.dot(Cinvbi,vectorbi))+np.log(np.linalg.det(Covbi))
                        chi2nomar = np.dot(diff,step1)
                    
                        chi2 = chi2mar+chi2nomar
                        
                #print(k_junc_high,ntry)

                if np.isnan(chi2):

                    chi2 = chi2_original
                    #print('chi2nan = ' + str(chi2))  
                    #print('theta = ')
                    #print(theta)
                    #print(time.time()-t0)
                    #print('ntry = '+ str(ntry))  
            else:
                    chi2 = chi2_original

        if withPlanck:
            return -0.5* ( chi2 + (rd-rs(Om,h,f_fid))**2/sigma_rd**2 )
        else:
        	return -0.5* chi2


def lnprob(theta, xdata, ydata, Cinv, free_para, fix_para,bounds,Om_fid,f_fid,binning=False,TableNkmu=None, window=True,dataQ=None,withBisp=False,masktriangle=None,Bispdata=None,withPlanck=False):
   
    """ Computes the log of the probability (logprior + loglike)
    Inputs
    ------
    theta : array of all parameters
    xdata : the k on which the data is evaluated (concatenated for multipoles)
    ydata : the data concatenated for multipoles
    Cinv : the inverse of the covariance
    free_para : list of boolean corresponding to the parameters to be varied
    fix_para : list of the values for the fixed parameters
    fiducial : the value of cosmological parameters on the fiducial (for AP effect)
    bounds : 2d array with the [min,max] values for the parameters
    interpolation_grid : the interpolation of the power spectra on the grid
    binning : whether to take into account the discreteness of the data power spectra. Must provide a number of modes per bin in TableNkmu
    TableNkmu : the number of modes per (k,mu) bin. Expected to be of the type kmean, mucentral, Nkmu
    window : whether to take into account the window function (when not a periodic square box, such as Lightcone). Must provide the data for the Q matrices (see e.g. section 4.3.2 of arXiv:1706.02362)
    dataQ : the Q matrices (multipoles of the window function) taken from Florian. Of the form stab,dummy,Q0,Q2,Q4,Q6,Q8,mueff
    Outputs
    ------
    the value of (logprior + loglike)
    
    """    
    lp  =  lnprior(theta, free_para, fix_para,bounds)
    
    if np.isfinite(lp) == False :
        dummy  =  -np.inf
        
    dummy  =  lp + lnlike(theta, xdata, ydata, Cinv, free_para, fix_para,bounds,Om_fid, f_fid,binning=binning,TableNkmu=TableNkmu, window=window,dataQ=dataQ,withBisp=withBisp,masktriangle=masktriangle,Bispdata=Bispdata,withPlanck=withPlanck)

    return dummy



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

    defaultgrid = sys.argv[4]
    
    series_cosmo = dfcosmo.loc[simtype]
    if defaultgrid == 1:
        gridname = series_cosmo.log['gridname']	 
    else:
	try:
	    sys.argv[5]
	except:
	    raise Exception("You asked for a non-default grid but didn't give its name as an argument!")
	gridname = sys.argv[5] 
    
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

    withPlanck = False 


    marg_gaussian = False
    
    withBisp = False 
 

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

    if ZONE != '':    
        dataQ = np.loadtxt(opa.join(INPATH,'Window_functions/dataQ_%s.txt'%ZONE)).T 
    elif 'ChallengeQuarter' in simtype:
        dataQ = np.loadtxt(opa.join(INPATH,'Window_functions/dataQ_ChallengeQuarter.dat')).T
    if 'Challenge' in simtype:
        Full_Cov = np.loadtxt(opa.join(INPATH,'Covariance/Cov%s%s.dat'%(simtype,ZONE)))
    else: 
        Full_Cov = np.loadtxt(opa.join(INPATH,'Covariance/Cov%s%sdata.dat'%(simtype,ZONE)))
    
    print("here")    
    #bra
    runtype = simtype+ZONE
    
    Bispdata = [] 
    masktriangle = []
    if withBisp:
        runtype += 'withBispkmax%s'%kmaxbisp
        Full_Cov = np.loadtxt(opa.join(INPATH,'Covariance/Cov%s%s_Bisp.dat'%(simtype,ZONE)))   
        Q1,Q2,Q3,Bispdata = np.loadtxt(opa.join(INPATH,'DataSims/Bispred_LightConeHector_%s_%s.dat'%(ZONE,boxnumber))).T
        masktriangle = (Q1>=kminbisp)&(Q1<=kmaxbisp)&(Q1<=Q2+Q3)&(Q1>=abs(Q2-Q3))&(Q2>=kminbisp)&(Q2<=kmaxbisp)&(Q3>=kminbisp)&(Q3<=kmaxbisp)
    #print(masktriangle.shape)
    window = True
    binning = False
    TableNkmu = None
    
    

    lnAsmin,lnAsmax,Ommin,Ommax,hmin,hmax,interpolation_grid = get_grid(gridname,nbinsAs=100,withBisp=withBisp)    
    print("got grid!")    
##############################
###  Priors ###################
#############################

    #### The uniform prior on the b_i#####
    bmin = -200
    bmax = 200

    # We require b_1>0
    bmintab = [0] +[bmin]*10
    bmaxtab = [bmax]*11
    

    ##### The bounds for the minimizer ######
    bfmaxtab = np.concatenate([[lnAsmax,Ommax,hmax],bmaxtab])
    bfmintab = np.concatenate([[lnAsmin,Ommin,hmin],bmintab])

    bounds = zip(bfmintab,bfmaxtab)

    ##### Initial guess for the b_i #####
    inipos = np.array([1.85 ,  -2.62623719,  -0.39661384,   4.21514113,
         8.36786486, -29.68630616,   1.03528956, 0,
        40.00717862,  0,   100])

    ##### Guess for the \sigma, to help with the initial position of walkers #####
    onesigma = np.array([  2.48736140e-01,   4.40317511e-02,   2.65820186e-02,
          3.13222972e-01,   1.13467462e+01,   3.17680821e+01,
          1.02673261e+01,   5.73925336e+01,   3.33253123e+01,
          2.49030346e+01,   5.73415031e+01,   1.38582379e+02,
          1.40667700e+02,1.40667700e+02]) 

    
    kmaxtab = [KMAX]

    print("Starting process")
    sys.stdout.flush()
    kmaxname = ['kmax%s'%kmax for kmax in kmaxtab]

    for boxnumber in [boxnumber]:
	if 'ChallengeQuarter' in simtype:
		simtype = simtype[:-2]
        kPS,PSdata,_ = np.loadtxt(opa.join(INPATH,'DataSims/ps1D_%s%s_%s.dat'%(simtype,ZONE,boxnumber))).T

        for indexk,kmax in enumerate(kmaxtab):    

            xdata = kPS[(kPS<kmax)&(kPS>kmin)]
            ydata = PSdata[(kPS<kmax)&(kPS>kmin)]
            indexkred = np.argwhere((kPS<kmax)&(kPS>kmin))[:,0]
            if withBisp:
                indextriangle = np.argwhere(masktriangle)[:,0]+kPS.shape[0]
                indexkred = np.concatenate([indexkred,indextriangle])
            print("triangle length is", np.argwhere(masktriangle)[:,0].shape)        
            Covred = Full_Cov[indexkred[:,None],indexkred]

            Cinv = np.linalg.inv(Covred)
    
    #################################
    ## Setting up the fit ###########
    #################################    
    
            all_true  =  np.concatenate(([lnAs_fid, Om_fid, h_fid],[1.9]+10*[0]))
            all_name  =  np.concatenate(([r'$A_s$',r'$\Omega_m$',r'$h$'],[r'$b_%s$'%i for i in range(len(inipos))]))
            free_para  =  [True,True,True,True,True,False,True,False,False,False,False,False,False,withBisp]
            
            nparam = len(free_para)
            
            
            # if free_para is false read the value in fix_para
            fix_para  =  all_true

            ndim  =  sum(free_para)
            
            free_true = all_true[free_para]
            free_name = all_name[free_para]



    #################################
    ## Find maximum likelihood ######
    #################################

        t0 = time.time()    
        chi2  =  lambda theta: -2 * lnlike(theta,xdata, ydata, Cinv, free_para, fix_para,bounds,Om_fid,f_fid,marg_gaussian=marg_gaussian,binning=binning,window=window,withBisp=withBisp,dataQ=dataQ,TableNkmu=TableNkmu,Bispdata=Bispdata,masktriangle=masktriangle,withPlanck=withPlanck)
        print("computed minimum!", time.time() - t0)
        sys.stdout.flush() 

        result  =  op.minimize(chi2, all_true,method = 'SLSQP',bounds = bounds,options = {'maxiter':100})

        all_ml  =  result["x"]
        
        free_ml = all_ml[free_para]

        minchi2  =  result["fun"]
        
        print(free_ml)
        print("ap params", get_AP_param(free_ml[1], Om_fid))
        if type(masktriangle) == type(None):
            dof = len(xdata) - ndim
        else:
            dof = len(xdata) + sum(masktriangle) - ndim
        print('minchi2 = ' + str(minchi2), dof)
        np.savetxt(opa.join(OUTPATH,"minchi2%sbox_%skmax_%s.txt")%(runtype,boxnumber,kmax),np.concatenate([free_ml,[minchi2,dof]]))
    
    ###################################
    ## run MCMC #######################
    ###################################

       # Set up the sampler.
 

    Nchains  =  4
    nwalkers  =  4*nparam
    fidpos = np.concatenate([ [ lnAs_fid,   Om_fid,   h_fid],  free_ml[3:]])


    # Start MCMC
    t0 = time.time()
    temperature  =  1.
    minlength  =  4000
    ichaincheck  =  50
    ithin  =  1
    epsilon  =  0.03
    # Set up the sampler.
    pos = []
    sampler = []
    rstate  =  np.random.get_state()

    print("ndim  =  ", ndim)
    print("Nchains  =  ", Nchains)

    for jj in range(0, Nchains):

        initialpos = []
        for ii in xrange(nwalkers):
            accepted  =  False
            while (not accepted):
                trialfiducial  =  np.random.normal(loc = free_ml,scale =  temperature*onesigma[free_para])
                accepted  =  np.isfinite(lnprior(trialfiducial, free_para, fix_para,bounds))
            if accepted:
                initialpos.append(trialfiducial)
        pos.append(initialpos)
        sampler.append(emcee.EnsembleSampler(nwalkers, ndim, lnprob,a = 1.15, args = (xdata, ydata, Cinv, free_para, fix_para,bounds,Om_fid,f_fid),kwargs={'marg_gaussian':marg_gaussian,'binning':binning,'window':window,'withBisp':withBisp,'dataQ':dataQ,'masktriangle':masktriangle,'TableNkmu':TableNkmu,'Bispdata':Bispdata,'withPlanck':withPlanck},threads = 4))
        
    np.save(opa.join(OUTPATH,"inipos%sbox_%skmax_%s")%(runtype,boxnumber,kmax),np.array(pos))
    # Start MCMC
    print("Running MCMC...")

    withinchainvar  =  np.zeros((Nchains,ndim))
    meanchain  =  np.zeros((Nchains,ndim))
    scalereduction  =  np.arange(ndim,dtype = np.float)
    for jj in range(0, ndim):
        scalereduction[jj]  =  2.
    sys.stdout.flush()
    itercounter  =  0
    chainstep  =  minlength
    loopcriteria  =  1
    t1 = time.time()
    for jj in range(0, Nchains):
        try:
            progress = np.load(opa.join(OUTPATH,"ChainsMidway/samplerchainmid%sbox_%skmax_%srun_%s")%(runtype,boxnumber,kmax,jj))
            print("found old chains!, size is ", progress.shape)
            sampler[jj].chain = progress
        except:
            pass

    while loopcriteria:
        
        itercounter  =  itercounter + chainstep
        print("chain length  = ",itercounter," minlength  = ",minlength)
        samplesJG = []
        for jj in range(0, Nchains):
            # Since we write the chain to a file we could put storechain = False, but in that case
            # the function sampler.get_autocorr_time() below will give an error
            c = 0
	    for result in sampler[jj].sample(pos[jj], iterations = chainstep, rstate0 = rstate, storechain = True, thin = ithin):
                pos[jj]  =  result[0]
                chainchi2  =  -2.*result[1]
                rstate  =  result[2]
         #       if c%10==0:
		    #print(time.time() - t0, c)
	#	c+=1
	    t0 = time.time()
    
            # we do the convergence test on the second half of the current chain (itercounter/2)
            chainsamples  =  sampler[jj].chain[:, itercounter/2:, :].reshape((-1, ndim))
            #print("len chain  =  ", chainsamples.shape)
            withinchainvar[jj]  =  np.var(chainsamples, axis = 0)
            meanchain[jj]  =  np.mean(chainsamples, axis = 0)
            samplesJG.append(chainsamples)
            np.save(opa.join(OUTPATH,"ChainsMidway/samplerchainmid%sbox_%skmax_%srun_%s")%(runtype,boxnumber,kmax,jj),sampler[jj].chain[:,::10,:]) 
            np.save(opa.join(OUTPATH,"ChainsMidway/lnlikechainmid%sbox_%skmax_%srun_%s")%(runtype,boxnumber,kmax,jj),sampler[jj].lnprobability[:20,::10])
        print(time.time() - t1)
#	sys.stdout.flush()
	scalereduction  =  gelman_rubin_convergence(withinchainvar, meanchain, itercounter/2, Nchains, ndim)
        print("scalereduction  =  ", scalereduction)
        sys.stdout.flush() 
        loopcriteria  =  0
        for jj in range(0, ndim):
            if np.absolute(1-scalereduction[jj]) > epsilon:
                loopcriteria  =  1

        chainstep  =  ichaincheck

    print("Done.")
    trun = time.time()-t0
    print(trun)

    

    # Print out the mean acceptance fraction. In general, acceptance_fraction
    # has an entry for each walker so, in this case, it is a 250-dimensional
    # vector.
    for jj in range(0, Nchains):
        print("Mean acceptance fraction for chain ", jj,": ", np.mean(sampler[jj].acceptance_fraction))
        np.savetxt(opa.join(OUTPATH,'AcceptanceFr%sbox_%skmax_%srun_%s.dat'%(runtype,boxnumber,kmax,jj)),[np.mean(sampler[jj].acceptance_fraction)])
# Estimate the integrated autocorrelation time for the time series in each
    # parameter.
    #print("Autocorrelation time:", sampler.get_autocorr_time())
    
    burnin  =  1000
    
    for jj in range(Nchains):
        np.save(opa.join(OUTPATH,"samplerchain%sbox_%skmax_%srun_%s")%(runtype,boxnumber,kmax,jj),sampler[jj].chain[:,::5,:])
        np.save(opa.join(OUTPATH,"lnlikechain%sbox_%skmax_%srun_%s")%(runtype,boxnumber,kmax,jj),sampler[jj].lnprobability[:,::5])


    ###################################
    ## Compute the quantiles ##########
    ###################################

    mcmc_array  =  map(lambda v: (v[1], v[2]-v[1], v[1]-v[0], v[4]-v[1], v[1]-v[3]), zip(*np.percentile(np.array(samplesJG).reshape((-1,ndim)), [15.86555, 50, 84.13445, 2.2775, 97.7225], axis = 0)))


    np.savetxt(opa.join(OUTPATH,"mcmcarray%sbox_%skmax_%s.txt")%(runtype,boxnumber,kmax),mcmc_array)
