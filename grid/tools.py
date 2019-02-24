import os
import time
import hashlib
import configobj as cfg
import scipy
import numpy as np
import sys
import subprocess
import logging
from glob import glob

THIS_PATH = os.path.dirname(os.path.abspath(__file__))

# Expects the EFT code compiled here!!!
EFT_PATH = os.path.abspath(os.path.join(THIS_PATH, 'RedshiftBiasEFTwithFFT'))

metafilpath = os.path.join(EFT_PATH, 'metafil')

# Gets the metafil, some module with file utilities, and appends it to PATH
if os.path.isdir(metafilpath):
    if metafilpath not in sys.path:
        sys.path.append(metafilpath)
    import metafil
else:
    raise Exception('Module not found at ' + metafilpath)

# zbEFT file naming convention
FBSNM = 'PowerSpectra'
FEFT = '1loop'
FLIN = 'Linear'
Trianglepath = os.path.abspath(
    os.path.join(EFT_PATH, 'TrianglesConfiguration'))

# Parameters for config files.
# We divide the keys in 3 parts
# This is for CLASS.
# We do not use prior from BBN in order to explore the full parameter space (h, Om, As)
# (with Omb/Omc fixed). The relative difference in the matter power spectrum is insignificant anyway...

DEFAULTCONFIG_CLASS = ['ln10^{10}A_s', 'n_s', 'h', 'omega_b', 'omega_cdm',
                       'output', 'P_k_max_h/Mpc', 'root', 'headers', 'format', 'YHe', 'z_pk']

# This config is to for the EFT grid
DEFAULTCONFIG_zbEFT = ['knl', 'km', 'nbar',
                       'PathToLinearPowerSpectrum', 'PathToFolderOutput', 'PathToFolderRD',
                       'PathToFolderCosmoRef', 'ComputePowerSpectrum', 'UseCosmoRef',
                       'ImportResummationMatrix', 'ExportResummationMatrix', 'ComputeBispectrum',
                       'EpsRel_IntegrBispectrumAP', 'PathToTriangles',
                       'z_pk', 'omega_b', 'omega_cdm', 'ln10^{10}A_s', 'n_s', 'h',
                       'EpsAbs_NoCosmoRef', 'EpsRel_NoCosmoRef', 'EpsAbs_YesCosmoRef',
                       'EpsRel_YesCosmoRef', 'aperp', 'apar']
# For running with a window function
DEFAULTCONFIG_zbEFTw = ['outpath', 'basename', 'pid', 'CLASS_configf',
                        'zbEFT_configf', 'zbEFTw_configf', 'logfile',
                        'CLASS_path', 'CLASS_exe', 'CLASS_pre', 'zbEFT_path', 'zbEFT_exe',
                        'DM', 'kren']


def cH(Om, a):
    """Dimensionless Hubble in conformal time"""
    return np.sqrt(Om / a + a**2 * (1 - Om))


def DgN(Om, a):
    """Linear growth factor"""
    return 5. / 2 * Om * cH(Om, a) / a * scipy.integrate.quad(lambda x: cH(Om, x)**-3, 0, a)[0]


def fN(Om, a):
    """d ln D / d ln a"""
    return (Om * (5 * a - 3 * DgN(Om, a))) / (2. * (a**3 * (1 - Om) + Om) * DgN(Om, a))


def Classtomultipoles(setk, Pkclass, Om, zpk):
    """This are the linear (Kaiser) multipoles if b_1 = 1"""
    kmulti = np.concatenate([setk, setk, setk])
    f = fN(Om, 1. / (+zpk))
    P0k = (1 + (2 * f) / 3. + f**2 / 5.) * Pkclass
    P2k = ((4 * f * (7 + 3 * f)) / 21.) * Pkclass
    P4k = ((8 * f**2) / 35.) * Pkclass
    return np.array([kmulti, np.concatenate([P0k, P2k, P4k])])


def multipolestoClass(kmulti, Pmulti, Om, zpk):
    """This is pretty useles"""
    idxonethird = int(len(kmulti) / 3)
    kclass = kmulti[:idxonethird]
    P0k = Pmulti[:idxonethird]
    # P2k = Pmulti[idxonethird:2 * idxonethird]
    # P4k = Pmulti[2 * idxonethird:3 * idxonethird]
    f = fN(Om, 1. / zpk)
    Pkclass = P0k / (1 + (2 * f) / 3. + f**2 / 5.)
    return np.array([kclass, Pkclass])


def get_timestamp(prefix='', suffix=''):
    t = time.mktime((2019, 1, 1, 0, 0, 0, 0, 0, 0))
    if isinstance(prefix, (list, np.ndarray)):
        prefix = '_'.join(prefix)
    if isinstance(suffix, (list, np.ndarray)):
        suffix = '_'.join(suffix)
    return "%s_%s_%s" % (str(prefix), str(round(time.time()-t, 4)), str(suffix))


def make_hash(config_dict):
    hasher = hashlib.sha1()
    hasher.update(str(config_dict).encode("utf-8"))
    return hasher.hexdigest()[:7]


def combine_hashes(configs):
    return '-'.join([str(make_hash(c)) for c in configs])


def setup_outputs(config_class, config_zbEFT, config_zbEFTw):
    """Create the name for the unique output directory using this structure:
    path/basename_timestamp_githash_classconfighash_zbeftconfighash_zbeftwconfighash/ """
    git = metafil.GitEnv()
    basestr = get_timestamp(config_zbEFTw['basename'],
                            [git.get_hash(7),
                             combine_hashes([config_class, config_zbEFT, config_zbEFTw])])
    if config_zbEFTw['pid']:
        basestr = basestr + '_' + str(config_zbEFTw['pid'])
    if basestr in config_zbEFTw['outpath']:
        raise IOError('It seems like you are reusing an ini file from a previous run. ' +
                      'Please create a copy and set the outpath parameter to a parent folder!')
    config_zbEFTw['outpath'] = os.path.abspath(
        os.path.join(config_zbEFTw['outpath'], basestr))

    # Replace the paths in the config so that they can be easily accessed as absolute
    # Check what has been provided - prepend folder if none is present
    [
        config_zbEFTw['CLASS_exe'],  # Where to find the class executable
        config_zbEFTw['CLASS_pre'],
        config_class['root'],  # Where to save the class linear pk
        config_zbEFTw['zbEFT_exe'],  # Where to find the zbeft executable
        # General output directory for zbeft
        config_zbEFT['PathToFolderOutput'],
        config_zbEFTw['logfile'],  # What to save the wzbeft log file to
        # Where to find/save intermediate files: resummation data
        config_zbEFT['PathToFolderRD'],
        # Where to find/save intermediate files: resummation data
        config_zbEFT['PathToFolderCosmoRef']
    ] = safe_prepend_folder([
        config_zbEFTw['CLASS_path'],
        config_zbEFTw['CLASS_path'],
        config_zbEFTw['outpath'],
        config_zbEFTw['zbEFT_path'],
        config_zbEFTw['outpath'],
        config_zbEFTw['outpath'],
        config_zbEFTw['outpath'],
        config_zbEFTw['outpath']
    ], [
        config_zbEFTw['CLASS_exe'],
        config_zbEFTw['CLASS_pre'],
        config_class['root'],
        config_zbEFTw['zbEFT_exe'],
        '',  # PathToFolderOutput should always be set automatically!!
        config_zbEFTw['logfile'],
        config_zbEFT['PathToFolderRD'],
        config_zbEFT['PathToFolderCosmoRef']
    ])

    # Finish the CLASS power spectrum file name
    # hardcoded, because this is what class does!
    # Finish the CLASS power spectrum file name
    config_zbEFT['PathToLinearPowerSpectrum'] = config_class['root'] + \
        'pk.dat'  # hardcoded, because this is what class does!

    # Make the output folder if it doesn't yet exist
    try:
        if not os.path.isdir(config_zbEFTw['outpath']):
            os.makedirs(config_zbEFTw['outpath'])
        else:
            os.makedirs(config_zbEFT['PathToFolderOutput'])
    except IOError:
        print('Cannot create directory: %s' % config_zbEFTw['outpath'])
    else:
        logging.info(
            'Created new output directory and subdirectories: %s' % config_zbEFTw['outpath'])

    # If needed create folders for resummation data
    if config_zbEFT['ExportResummationMatrix'] == 'yes':
        if os.path.isdir(config_zbEFT['PathToFolderRD']):
            print("Resummation data folder already exists. Any files will be overwritten.\n %s"
                  % config_zbEFT['PathToFolderRD'])
        else:
            os.makedirs(config_zbEFT['PathToFolderRD'])

    return basestr


def make_configfiles(configs,
                     filekeys=['CLASS_configf',
                               'zbEFT_configf', 'zbEFTw_configf'],
                     indfmeta=2):
    """Writes the config files as output... This is pretty terrible Python.
    And wtf is indfmeta?"""
    for i, c in enumerate(configs):
        found = True
        for fk in filekeys:
            if fk not in c.keys():
                found = False
        if found:
            break
    if not found:
        raise Exception(
            'Configuration for the zbEFT wrapper not found in input!')
    if filekeys[i] not in configs[i].keys():
        raise Exception('The order of the inputs must match!')

    for c, f in zip(configs, filekeys):
        c.filename = os.path.join(
            configs[indfmeta]['outpath'], configs[indfmeta][f])
        c.write()
        [configs[indfmeta][f]] = safe_prepend_folder(
            configs[indfmeta]['outpath'], configs[indfmeta][f])

    return


def get_config(bigconfig_file=None, bigconfig=None, cat=False):
    # Checked and modified
    """ Splits the configuration into 3: one for CLASS, one for the EFT code,
        one for the wrapper. Whatever is unset should be set to the default value.
        Inputs
        ------
        EITHER
        bigconfig_file : str
            path to the full config file
        OR
        bigconfig : dict or cfg.ConfigObj
            the full configuration in a dict or ConfigObj
        cat: Tells us if we want to concatenate the output.
    """

    if bigconfig_file is not None:
        if bigconfig is not None:
            raise Exception('get_config only takes one input!')
        if not os.path.isfile(os.path.abspath(bigconfig_file)):
            raise IOError(bigconfig_file + ' file not found!')
        bigconfig = cfg.ConfigObj(bigconfig_file)
        # bigconfig.filename = bigconfig_file

    saved = []

    cclass = cfg.ConfigObj()
    for key in DEFAULTCONFIG_CLASS:
        try:
            cclass[key] = bigconfig[key]
        except KeyError:
            print("Please set %s in the config file!" % key)
        saved.append(key)

    czbEFT = cfg.ConfigObj()
    for key in DEFAULTCONFIG_zbEFT:
        try:
            czbEFT[key] = bigconfig[key]
        except KeyError:
            print("Please set %s in the config file!" % key)
        saved.append(key)

    czbEFTw = cfg.ConfigObj()
    for key in DEFAULTCONFIG_zbEFTw:
        try:
            czbEFTw[key] = bigconfig[key]
        except KeyError:
            print("Please set %s in the config file!" % key)
        saved.append(key)

    czbEFTw['pid'] = np.random.randint(0, 500)

    if not cat:
        return cclass, czbEFT, czbEFTw
    else:
        # return cclass.merge(czbEFT).merge(czbEFTw)
        cclass.merge(czbEFT)
        cclass.merge(czbEFTw)
        return cclass


def prepend_folder(paths, fnames):
    if not isinstance(fnames, list):
        fnames = [fnames]
    if not isinstance(paths, list):
        paths = [paths] * len(fnames)
    elif len(paths) != len(fnames):
        raise Exception("Can't combine path and file lists of lengths %d and %d!" %
                        (len(paths), len(fnames)))
    for i, [p, fn] in enumerate(zip(paths, fnames)):
        fnames[i] = os.path.abspath(os.path.join(p, fn))
    return fnames


def safe_prepend_folder(paths, fnames):
    for i, fn in enumerate(fnames):
        if os.path.isdir(os.path.dirname(fn.strip('/'))):
            paths[i] = ''
    return prepend_folder(paths, fnames)


def runcommand(command, logfile, outfile=None, cwd='.'):
    """Run command, with its arguments, in the cwd directory.
    Prints the output to logfile, and creates outfile which is checked"""
    with open(logfile, "wb") as out:
        process = subprocess.Popen(command, stdout=out, stderr=out, cwd=cwd)
        try:
            process.wait()
        except KeyboardInterrupt as e:
            process.kill()
            raise e

    if outfile is not None:
        if len(glob(outfile)) == 0:
            errmsg = "Command: %s failed for unknown reasons. \
            The following expected file not found: %s." % (' '.join(command), outfile)
            raise Exception(errmsg)
    return


def read_file(filenamebase, multipole, path='./', verb=0):

    filename = os.path.abspath(
        os.path.join(path, filenamebase + '_l' + str(multipole) + '.dat'))
    if verb:
        print('reading from ' + filename)

    # Get the column names first
    with open(filename, 'r') as f:
        cols = f.readline().split()

    # Check that all is as expected
    if '#' not in cols:
        raise IOError(
            'A header beginning with # is expected in file: ' + filename)
    elif cols.index('k') != 1:
        raise IOError('Expecting k as first column in file: ' + filename)
    cols.pop(0)
    cols.pop(0)

    # Get the data
    data = np.loadtxt(filename)

    return data[:, 0], data[:, 1:], cols


def into_1arr(ls, ks, mlps, mlps_lin, nms, nms_lin):

    # Precalculate lengths
    nks = [len(k) for k in ks]
    nktot = sum(nks)
    nncol = None
    nlcol = None
    for ml, m in zip(mlps_lin, mlps):
        if nncol is None or nlcol is None:
            nlcol = len(ml[0, :])
            nncol = len(m[0, :])
        elif nlcol + nncol != len(ml[0, :]) + len(m[0, :]):
            raise IOError(
                "Multipoles have different number of columns. That won't work!")
        elif nlcol != len(nms_lin) or nncol != len(nms):
            raise IOError(
                "Term array sizes differ from lists of columns. That won't work!")

    # Make the empty arrays in advance
    mpls = np.zeros(nktot)
    kvals = np.zeros(nktot)
    terms = np.zeros((nktot, nlcol + nncol))

    for nk, l, k, mpl, mpl_lin in zip(nks, ls, ks, mlps, mlps_lin):
        mpls[nk * l / 2:nk * (l / 2 + 1)] = l
        kvals[nk * l / 2:nk * (l / 2 + 1)] = k
        terms[nk * l / 2:nk * (l / 2 + 1), :nlcol] = mpl_lin
        terms[nk * l / 2:nk * (l / 2 + 1), nlcol:] = mpl

    return mpls, kvals, terms, nms_lin + nms

# Check the multipoles - not foolproof, but at least something


# def check_mlps(mlps):
#     """This is pretty useless, comment it out"""
#     for l in mlps:
#         if l not in MULTIPOLES.keys():
#             raise Exception('We can only take multipoles of ' +
#                             str(MULTIPOLES.keys()) + '!')
#     return True

# See how many times the kvals reverse -> this + 1 should tell you how many multipoles are stored in a vector


def count_mlps(kvals):
    return np.sum(np.diff(kvals) < 0) + 1

# Extract only one of the repeated k-arrays (return the index - same assumption as above)


def get_k_once(kvals, index=False):
    ends = [len(kvals)]
    ends[:0] = [i + 1 for i in np.nonzero(np.diff(kvals) < 0)[0]]
    index_once = np.arange(ends[0])
    if index:
        return index_once
    else:
        return kvals[index_once]


def read_zbEFT(basepath, mlps=[0, 2, 4], eft=FEFT, lin=FLIN):

    terms = [None] * len(mlps)
    terms_lin = [None] * len(mlps)
    kvals_nl = [None] * len(mlps)
    kvals_lin = [None] * len(mlps)

    names = None
    nms_l = None
    for i, mlp in enumerate(mlps):

        kvals_lin[i], terms_lin[i], names_lin = read_file(basepath + lin, mlp)
        kvals_nl[i], terms[i], names_nl = read_file(basepath + eft, mlp)

        if names is None or nms_l is None:
            names = names_nl
            nms_l = names_lin
        elif np.any(kvals_nl[i] != kvals_lin[i]) or np.any(nms_l != names_lin) or np.any(names != names_nl):
            raise IOError(
                'Mismatch in the k or column names in the files at' + basepath + '!')

    return into_1arr(mlps, kvals_lin, terms, terms_lin, names_nl, nms_l)
