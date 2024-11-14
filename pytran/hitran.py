

from __future__ import print_function, division


from pytran.hitran_utils import qtips, get_molecule_nisops, get_iso_mass

import numpy as np


__all__ = ['lorentzian_profile', 'doppler_profile', 'voigt_profile', 'read_hitran2012_parfile',
               'calculate_hitran_xsec']


def lorentzian_profile(wavn, S, gamma, wavn0):
    '''
    Calculate a Lorentzian absorption profile.

    Parameters
    ----------
    wavn : ndarray
        The array of wavenumbers at which to sample the profile.
    S : float
        The absorption line "strength" (the integral of the entire line is equal to S).
    gamma : float
        The linewidth parameter of the profile.
    wavn0 : float
        The center position of the profile (in wavenumbers).

    Returns
    -------
    L : ndarray
        The sampled absorption profile.
    '''
    L = (S / np.pi) * gamma / ((wavn - wavn0)**2 + gamma**2)
    return(L)
    

def doppler_profile(wavn, S, alpha, wavn0):
    '''
    doppler profile
    '''
    G = (S / (np.sqrt(np.pi)*alpha)) * np.exp(-((wavn-wavn0) / alpha)**2)
    return(G)
    

def voigt_profile(wavn, S, alpha, gamma, wavn0):
    '''
    voigt profile
    '''
    from scipy.special import wofz
    z = ((wavn-wavn0) + 1j*gamma)/(alpha)
    V = S/(np.sqrt(np.pi)*alpha)*wofz(z).real
    return(V)


def read_hitran2012_parfile(filename, wavemin=None, wavemax=None, Smin=None, isolist=None):
    """
    Given a HITRAN2012-format text file, read in the parameters of the molecular absorption features.

    Parameters
    ----------
    filename : str
        The filename to read in.

    Return
    ------
    linelist : dict
        The dictionary of HITRAN linelist for the molecule.
    """

    if filename.endswith('.zip'):
        import zipfile
        import os
        zipf = zipfile.ZipFile(filename, 'r')
        (object_name, ext) = os.path.splitext(os.path.basename(filename))
        print(object_name, ext)
        filehandle = zipf.read(object_name).splitlines()
    else:
        filehandle = open(filename, 'r')

    linelist = {
                'M':[],               ## molecule identification number
                'I':[],               ## isotope number
                'vc':[],              ## line center wavenumber (in cm^{-1})
                'S':[],               ## line strength, in cm^{-1} / (molecule m^{-2})
                'Acoeff':[],          ## Einstein A coefficient (in s^{-1})
                'gamma-air':[],       ## line HWHM for air-broadening
                'gamma-self':[],      ## line HWHM for self-emission-broadening
                'Epp':[],             ## energy of lower transition level (in cm^{-1})
                'N':[],               ## temperature-dependent exponent for "gamma-air"
                'delta':[],           ## air-pressure shift, in cm^{-1} / atm
                'Vp':[],              ## upper-state "global" quanta index
                'Vpp':[],             ## lower-state "global" quanta index
                'Qp':[],              ## upper-state "local" quanta index
                'Qpp':[],             ## lower-state "local" quanta index
                'Ierr':[],            ## uncertainty indices
                'Iref':[],            ## reference indices
                'flag':[],            ## flag
                'gp':[],              ## statistical weight of the upper state
                'gpp':[]              ## statistical weight of the lower state
                }

    #print('Reading "' + filename + '" ...')

    isochar_to_isonum = {'1':1,'2':2,'3':3,'4':4,'5':5,'6':6,'7':7,'8':8,'9':9,'0':10,'A':11,'B':12}

    if (isolist is not None):
        if isinstance(isolist,int):
            isolist = (isolist,)
    
    for line in filehandle:

        if (line[0] == '#'):
            continue

        if (len(line) < 160):
            raise ImportError('The imported file ("' + filename + '") does not appear to be a HITRAN2012-format data file.')


        M = np.uint(line[:2])
        I = isochar_to_isonum[line[2]]
        vc = np.float64(line[3:15])
        S = np.float64(line[15:25])


        if (wavemin is not None) and (vc < wavemin):
            continue
        if (Smin is not None) and (S < Smin):
            continue
        if (isolist is not None) and (I not in isolist):
            continue
        if (wavemax is not None) and (vc > wavemax):
            break  # don't need to continue to save


        # use line
        linelist['M'].append(np.uint(M))
        linelist['I'].append(np.uint(I))
        linelist['vc'].append(np.float64(line[3:15]))
        linelist['S'].append(np.float64(line[15:25]))
        linelist['Acoeff'].append(np.float64(line[25:35]))
        linelist['gamma-air'].append(np.float64(line[35:40]))
        linelist['gamma-self'].append(np.float64(line[40:45]))
        linelist['Epp'].append(np.float64(line[45:55]))
        linelist['N'].append(np.float64(line[55:59]))
        linelist['delta'].append(np.float64(line[59:67]))
        linelist['Vp'].append(line[67:82])
        linelist['Vpp'].append(line[82:97])
        linelist['Qp'].append(line[97:112])
        linelist['Qpp'].append(line[112:127])
        linelist['Ierr'].append(line[127:133])
        linelist['Iref'].append(line[133:145])
        linelist['flag'].append(line[145])
        linelist['gp'].append(np.float64(line[146:153]))
        linelist['gpp'].append(np.float64(line[153:160]))

    if filename.endswith('.zip'):
        zipf.close()
    else:
        filehandle.close()

    for key in linelist:
        linelist[key] = np.array(linelist[key])

    return linelist


def calculate_hitran_xsec(linelist, wavenumbers, M=None, T=296., P=101325., Pref=101325., qmix=0.0, wreach=25.0):
    """
    Given the HITRAN linelist (line centers and line strengths) for a molecule, digitize the result into a spectrum of
    absorption cross-section in units of cm^2.

    Parameters
    ----------
    linelist : dict of ndarrays
        The HITRAN data corresponding to a given molecule.
    wavenumber : ndarray of float
        The wavenumber to compute the cross section
    M : int, optional 
        hitran molecule number. If not given, then read from linelist
    T : float
        Temperature in Kelvin (Tref = 296. K)
    P : float
        Pressure in same units as Pref (default is Pref in Pa)
    Pref : float
        Reference pressure, use to change pressure units (Pref = 1 atm = 101325 Pa)
    qmix : float
        mixing ration of species in abosolute units (e.g. 0. to 1.; default is 0.)
    wreach : float
        extend in wavenumber to compute individual lines (default is 25 cm^-1)

    Returns
    -------
    xsec : array_like
        The mean absorption cross-section, in cm^2 per molecule, evaluated at the wavenumbers given by input `wavenumbers`.
    """

    #use constants from hitran for consistency, not from scipy.constants
    kb = 1.3806488e-16        # erg K−1
    hck = 1.4388
    amu =  1.660539066e-24   # g
    c = 2.99792458e10        # cm s−1
    #print(kb, c, amu, h, hck)

    # initialize xsec
    xsec = np.zeros_like(wavenumbers)

    # check M (M to be removed)
    if M == None:
        unique_M = np.unique(linelist['M'])
        if len(unique_M) != 1:
            raise Exception('pytran.calculate_hitran_crosssection expects only one molecule in linelist')
        M = unique_M[0]

    # check for errors in Epp (was an issue in previous versions of HITRAN)
    if any(linelist['Epp']<0.0):
        raise ValueError("Need to check for error flags in Epp")

    # scale line strengths to temperature

    #print(hck, np.max(linelist['S']))
    hckt = hck*(1./T-1./296.)
    qratio = [qtips(296.0,M,i)/qtips(T,M,i) for i in range(1,get_molecule_nisops(M)+1)]
    linestrengths = [S*qratio[int(I-1)]*np.exp(-hckt*Epp)*(1.0-np.exp(-hck*vc/T))/(1.0-np.exp(-hck*vc/296.0))
                     for (vc,S,I,Epp) in zip(linelist['vc'],linelist['S'],linelist['I'],linelist['Epp'])]

    # pre-calculate doppler widths and lorentz widths
    vcs = linelist['vc'] + linelist['delta']*(P/Pref)
    masses = [get_iso_mass(M,i)*amu for i in range(1,get_molecule_nisops(M)+1)]
    alphas = [np.sqrt(2.0*kb*T/masses[int(I-1)]) for I in linelist['I']] * vcs/c
    gammas = ((1.0-qmix)*linelist['gamma-air']+qmix*linelist['gamma-self']) * (P/Pref) * (296./T)**linelist['N']

    #print(P, Pref, T, 296.)
    #for il, linecenter, linestrength, alpha, gamma in zip(range(len(vcs)), vcs, linestrengths, alphas, gammas):
    #    print(linecenter, linestrength, alpha, gamma) 

    #
    for linecenter, linestrength, alpha, gamma in zip(vcs, linestrengths, alphas, gammas):
        #print(linecenter, linestrength, alpha/np.sqrt(2.), gamma) 
        #L = lorentzian_profile(wavenumbers, linestrength, gamma, linecenter)
        #xsec += L
        #D = doppler_profile(wavenumbers, linestrength, alpha, linecenter)
        #xsec += D
        #V = voigt_profile(wavenumbers, linestrength, alpha, gamma, linecenter)
        #xsec += V

        ## Note: the quantity sum(L * dk) should sum to "S"!
        #print(linestrength, V.sum()*(wavenumbers[1]-wavenumbers[0]),  D.sum()*(wavenumbers[1]-wavenumbers[0]),  L.sum()*(wavenumbers[1]-wavenumbers[0]))

        i1 = np.searchsorted(wavenumbers, linecenter-wreach)
        i2 = np.searchsorted(wavenumbers, linecenter+wreach)
        #print(linecenter, i1, i2, len(wavenumbers))
        V = voigt_profile(wavenumbers[i1:i2], linestrength, alpha, gamma, linecenter)
        xsec[i1:i2] += V

    return xsec
