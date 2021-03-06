from __future__ import print_function

import os
os.environ["TMPDIR"] = "/tmp"
import emcee
import numpy as np
import scipy.stats
import scipy.optimize as op
from scipy.interpolate import interp1d
import pandas as pd
import os.path as opa
import time
from scipy import stats

###########################################
###  Globals ##############################
###########################################

# PATH GLOBALS & tests to check they are as expected

# Save this path (expect data and other needed executables present in relative paths!!!)
THIS_PATH  =  opa.dirname(__file__)

import APpowerspectraNkmu
import getcorrelationfunctionwithCF
import WindowFFTlog

# Data paths
INPATH = opa.abspath(opa.join(THIS_PATH,'input')) 
OUTPATH = opa.abspath(opa.join(THIS_PATH,'output')) 



def Hubble(Om,h,z):
    return h*((Om)*(1+z)**3.+(1-Om))**0.5

def DA(Om,h,z):
    r = scipy.integrate.quad(lambda x:1./Hubble(Om,h,x), 0, z)[0]
    return r/(1+z)  

def get_AP_param(Om,h,fiducial):
    lnsAsfid,Omfid,hfid = fiducial
        
    qperp  =  DA(Om,h,z_pk)/DA(Omfid,hfid,z_pk)
    qpar  =  Hubble(Omfid,hfid,z_pk)/Hubble(Om,h,z_pk)
    
    return qperp,qpar

def check_if_multipoles_k_array(setk):
    return setk[len(setk)/3] == setk[0]    
            


def get_grid(gridname,nbinsAs=100,nbins = 50):
    
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

    Plininterp = scipy.interpolate.RegularGridInterpolator((lnAstab,Omtab,htab),TablePlin.reshape((nbinsAs,nbins,nbins,TablePlin.shape[-2],TablePlin.shape[-1])))
    Ploopinterp = scipy.interpolate.RegularGridInterpolator((lnAstab,Omtab,htab),TablePloop.reshape((nbinsAs,nbins,nbins,TablePloop.shape[-2],TablePloop.shape[-1])))
    interpolations = [Plininterp,Ploopinterp]

        
    return lnAsmin,lnAsmax,Ommin,Ommax,hmin,hmax,interpolations
    
    
    
def computePS(cvals,datalin,dataloop,setkin,setkout,removesqsig=True):
    
    """ Computes the power spectra given the b_i and the EFT power spectra

        Inputs
        ------
        bvals : The values for the b_i
        datalin : the linear power spectra from the EFT, with shape (multipoles, b_i, k)
        dataloop : the loop power spectra from the EFT, with shape (multipoles, b_i, k)
        setkin : the values of k from the input power spectra (must match datalin and dataloop)
        setkout : the values of k for the output power spectra
        removesqsig: whether or not to remove the piece of P22 which is constant at low k

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
    P11int =  interp1d(setkin, datalin[0,-1])
    sigsq = 2*scipy.integrate.quad(lambda x: x**2/(2*np.pi)**2*(P11int(x))**2,setkin.min(),setkin.max())[0]      
        
    P0 = interp1d(setkin,np.dot(cvals,data0)+datalin0[0]+b1*datalin0[1]+b1*b1*datalin0[2] - 2*removesqsig*(-b1 + b2 + b4)**2*sigsq)(setkout)
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
    
    
def prepare_data_single(boxnumber,multipole,rmin=30):
    rD,xiD = np.loadtxt(opa.join(INPATH,'CFdata/dr12v6c_rmu_cmass_NS_z43z75000%s.%s'%(boxnumber,multipole))).T[:,1:]

    rD1,rD2,CovD = np.loadtxt(opa.join(INPATH,'CFdata/dr12v6c_rmu_cmass_NS_z43z75_bin5.%s.cov.ave'%multipole)).T
    nr = np.unique(rD1).shape[0]
    CovDij = CovD.reshape((nr,nr))[1:,1:]
    rdata = rD[rD >= rmin]
    xidata = xiD[rD >= rmin]
    maskr = np.argwhere(rD>=rmin)[:,0] 

    Covred = CovDij[maskr[:,None],maskr]
    incovdata = np.linalg.inv(Covred)
    return rdata,xidata,incovdata
    
def prepare_data(boxnumber,rmin=30):
    
    xidata = []
    invcovdata = []
    

    for multipole in ['mono','quad','four']:
        result = prepare_data_single(boxnumber,multipole,rmin=rmin)
        rdata = result[0]
        xidata.append(result[1])
        invcovdata.append(result[2])
    
    return rdata,np.array(xidata),np.array(invcovdata)
    
def getchi2multi(rdata,xiD,r,xi,incov):
    
    xir = interp1d(r,xi,kind='cubic')(rdata)
    
    modelX = xir - xiD
    
    step1 = np.einsum('ij,j->i',incov,modelX)
    
    chi2 = np.dot(modelX,step1)
    return chi2

        
                
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


def lnlike(theta, rdata,xidata,invcovdata, free_para, fix_para,bounds,fiducial, interpolation_grid,binning=False,TableNkmu=None, window=False,dataQ=None,**kwargs):
    
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
        lnAs,Om,h,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11 = match_para(theta, free_para, fix_para)
        
    # Import the power spectra interpolators on the grid
        Plininterp,Ploopinterp = interpolation_grid 
        
        kfull = Ploopinterp((lnAs,Om,h))[:,0]
    
        if check_if_multipoles_k_array(kfull):
            kfull = kfull[:len(kfull)/3] 
        
        Ploop = np.swapaxes(Ploopinterp((lnAs,Om,h)).reshape(3,len(kfull),22),axis1 = 1,axis2 = 2)[:,1:,:]
        Plin = np.swapaxes(Plininterp((lnAs,Om,h)).reshape(3,len(kfull),4),axis1 = 1,axis2 = 2)[:,1:,:]       
      
     
        # Compute the PS       
        valueb = np.array([b1,b2,b3,b4,b5,b6,b7,b8,b9,b10])        
        Pmodel = computePS(valueb,Plin,Ploop,kfull,kfull)
    
        #The AP parameters
        qperp,qpar = get_AP_param(Om,h,fiducial)
    
        if not binning:
            Pmodel = APpowerspectraNkmu.changetoAPnobinning(Pmodel,kfull,kfull,qperp,qpar)
        else:
            if type(TableNkmu) == type(None):
                raise Exception('You want to account for binning but forgot to provide a TableNkmu (array of shape (3,n)) obtained from the sims/ Can be found in input/TableNkmu')
            else : Pmodel = APpowerspectraNkmu.changetoAPbinning(Pmodel,kfull,kfull,qperp,qpar,TableNkmu)
    
        
        Pkin = np.concatenate(Pmodel)
        setkin = kfull
    
        r,xi0,xi2,xi4 = getcorrelationfunctionwithCF.Ptoxi(Pkin,setkin,n=2**12,extrap=True,k_junc_low=setkin.min(),k_junc_high=0.9*setkin.max())
        
        ximodel = np.array([xi0,xi2,xi4])
        
        if window:
            if type(dataQ) ==  type(None):
                raise Exception('You want to account for window function but forgot to provide a dataQ (array of shape (8,n)) obtained from the sims. Can be found in input/dataQ')
            else: ximodel = WindowFFTlog.transform_CF(np.array([xi0,xi2,xi4]),r,dataQ)
        
                    
        
    
        chi20 = getchi2multi(rdata,xidata[0],r,ximodel[0],invcovdata[0])
        chi22 = getchi2multi(rdata,xidata[1],r,ximodel[1],invcovdata[1])
        chi24 = getchi2multi(rdata,xidata[2],r,ximodel[2],invcovdata[2])
        
        chi2tot = (chi20+chi22+chi24)
        if np.isnan(ximodel).any():
            return -200000
        else:
            return -0.5*chi2tot


def lnprob(theta, rdata,xidata,invcovdata, free_para, fix_para,bounds,fiducial, interpolation_grid,binning=False,TableNkmu=None, window=False,dataQ=None,**kwargs):
   
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
        
    dummy  =  lp + lnlike(theta, rdata,xidata,invcovdata, free_para, fix_para,bounds,fiducial, interpolation_grid,binning=binning,TableNkmu=TableNkmu, window=window,dataQ=dataQ,**kwargs)

    return dummy        


                                                                                                    ###########################################
                                                                                                    ###  Main program  ########################
                                                                                                    ###########################################


if __name__ ==  "__main__":

    # Table of cosmological parameters according to seems

    dfcosmo = pd.read_csv(opa.join(INPATH,'DataFrameCosmosims.csv'),index_col=0)
    simtype = "CFAlbert"
    
    
    
    # Load the row that we are interested in
    series_cosmo = dfcosmo.loc[simtype]
    
    
    gridname = series_cosmo.loc['gridname']
    

    # COSMOLOGICAL GLOBALS: fiducial model (should match input sim data!)
    Om_fid  =  series_cosmo.loc['Omega_m']
    lnAs_fid = series_cosmo.loc['lnAs']
    h_fid  =  series_cosmo.loc['h']
    z_pk = series_cosmo.loc['z_pk']
    
    fiducial = [lnAs_fid,Om_fid,h_fid]

    #### Choice for the data #####
    #For lightcone simulations, need to specify north or south for now (later, merge the two but I'm missing the covariance for SGC
    ZONE = 'NGC'
    
    if ZONE != '':    
        dataQ = np.loadtxt(opa.join(INPATH,'Window_functions/dataQ_%s.txt'%ZONE)).T 
    
    runtype = simtype+ZONE
    rmin = 30
    boxnumber = 1
    lnAsmin,lnAsmax,Ommin,Ommax,hmin,hmax,interpolation_grid = get_grid(gridname)
    
    binning = False
    window = False
    
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
    inipos = np.array([2.0310343 ,  -2.62623719,  -0.39661384,   4.21514113,
         8.36786486, -29.68630616,   1.03528956, -32.39092667,
        40.00717862,   4.61905778,   100])

    ##### Guess for the \sigma, to help with the initial position of walkers #####
    onesigma = np.array([  2.48736140e-01,   4.40317511e-02,   2.65820186e-02,
          3.13222972e-01,   1.13467462e+01,   3.17680821e+01,
          1.02673261e+01,   5.73925336e+01,   3.33253123e+01,
          2.49030346e+01,   5.73415031e+01,   1.38582379e+02,
          1.40667700e+02,1.40667700e+02]) 

    
    rmintab = [rmin]

    print("Starting process")
    

    for boxnumber in [boxnumber]:


        for indexk,rmin in enumerate(rmintab):    

            rdata,xidata,invcovdata = prepare_data(boxnumber,rmin)
    
    #################################
    ## Setting up the fit ###########
    #################################    
    
            all_true  =  np.concatenate(([lnAs_fid, Om_fid, h_fid],[1.9]+10*[0]))
            all_name  =  np.concatenate(([r'$A_s$',r'$\Omega_m$',r'$h$'],[r'$b_%s$'%i for i in range(len(inipos))]))
            free_para  =  [True,True,True,True,True,True,True,True,True,True,False,False,False,False]
            
            nparam = len(free_para)
            
            
            # if free_para is false read the value in fix_para
            fix_para  =  all_true
            # create an array of free parameters
            counter  =  0
            for i in range(nparam):
                if free_para[i]  ==  True:
                    counter +=  1

            ndim  =  counter
            
            free_true = []
            free_name = []
            free_ml  =  np.arange(counter,dtype = np.float)

            counter  =  0;
            for i in range(nparam):
                if free_para[i]  ==  True:
                    free_true.append(all_true[i])
                    free_name.append(all_name[i])
                    counter +=  1


    #################################
    ## Find maximum likelihood ######
    #################################

        
        chi2  =  lambda theta: -2 * lnlike(theta, rdata,xidata,invcovdata, free_para, fix_para,bounds,fiducial, interpolation_grid,binning=binning,window=window,dataQ=dataQ)
    

        result  =  op.minimize(chi2, all_true,method = 'SLSQP',bounds = bounds,options = {'maxiter':100})

        all_ml  =  result["x"]
        
        free_ml = all_ml[free_para]

        minchi2  =  result["fun"]
        
        dof = len(xidata.reshape(-1)) - ndim
    
        np.savetxt(opa.join(OUTPATH,"minchi2%sbox_%skmax_%s.txt")%(runtype,boxnumber,rmin),np.concatenate([free_ml,[minchi2,dof]]))

    ###################################
    ## run MCMC #######################
    ###################################

       # Set up the sampler.
 

    Nchains  =  4
    nwalkers  =  2*nparam
    fidpos = np.concatenate([ [ lnAs_fid,   Om_fid,   h_fid],  free_ml[3:]])


    # Start MCMC
    t0 = time.time()
    temperature  =  1.
    minlength  =  4000
    ichaincheck  =  50
    ithin  =  1
    epsilon  =  0.06
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
        sampler.append(emcee.EnsembleSampler(nwalkers, ndim, lnprob,a = 1.15, args = (rdata,xidata,invcovdata, free_para, fix_para,bounds,fiducial, interpolation_grid),kwargs={'binning':binning,'window':window,'dataQ':dataQ},threads = 1))
        
    np.save(opa.join(OUTPATH,"inipos%sbox_%skmax_%s")%(runtype,boxnumber,rmin),np.array(pos))
    # Start MCMC
    print("Running MCMC...")

    withinchainvar  =  np.zeros((Nchains,ndim))
    meanchain  =  np.zeros((Nchains,ndim))
    scalereduction  =  np.arange(ndim,dtype = np.float)
    for jj in range(0, ndim):
        scalereduction[jj]  =  2.

    itercounter  =  0
    chainstep  =  minlength
    loopcriteria  =  1
    while loopcriteria:
        
        itercounter  =  itercounter + chainstep
        print("chain length  = ",itercounter," minlength  = ",minlength)
        samplesJG = []
        for jj in range(0, Nchains):
            # Since we write the chain to a file we could put storechain = False, but in that case
            # the function sampler.get_autocorr_time() below will give an error
            for result in sampler[jj].sample(pos[jj], iterations = chainstep, rstate0 = rstate, storechain = True, thin = ithin):
                pos[jj]  =  result[0]
                chainchi2  =  -2.*result[1]
                rstate  =  result[2]
    
    
            # we do the convergence test on the second half of the current chain (itercounter/2)
            chainsamples  =  sampler[jj].chain[:, itercounter/2:, :].reshape((-1, ndim))
            #print("len chain  =  ", chainsamples.shape)
            withinchainvar[jj]  =  np.var(chainsamples, axis = 0)
            meanchain[jj]  =  np.mean(chainsamples, axis = 0)
            samplesJG.append(chainsamples)
            np.save(opa.join(OUTPATH,"ChainsMidway/samplerchainmid%sbox_%skmax_%srun_%s")%(runtype,boxnumber,rmin,jj),sampler[jj].chain[:20,::10,:])
        scalereduction  =  gelman_rubin_convergence(withinchainvar, meanchain, itercounter/2, Nchains, ndim)
        print("scalereduction  =  ", scalereduction)
        
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
        np.savetxt(opa.join(OUTPATH,'AcceptanceFr%sbox_%skmax_%srun_%s.dat'%(runtype,boxnumber,rmin,jj)),[np.mean(sampler[jj].acceptance_fraction)])
# Estimate the integrated autocorrelation time for the time series in each
    # parameter.
    #print("Autocorrelation time:", sampler.get_autocorr_time())
    
    burnin  =  1000
    
    for jj in range(Nchains):
        np.save(opa.join(OUTPATH,"samplerchain%sbox_%skmax_%srun_%s")%(runtype,boxnumber,rmin,jj),sampler[jj].chain[:,::5,:])


    ###################################
    ## Compute the quantiles ##########
    ###################################

    mcmc_array  =  map(lambda v: (v[1], v[2]-v[1], v[1]-v[0], v[4]-v[1], v[1]-v[3]), zip(*np.percentile(np.array(samplesJG).reshape((-1,ndim)), [15.86555, 50, 84.13445, 2.2775, 97.7225], axis = 0)))

    np.savetxt(opa.join(OUTPATH,"mcmcarray%sbox_%skmax_%s.txt")%(runtype,boxnumber,rmin),mcmc_array)


    
    



