

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


def read_hitran2012_parfile(filename, wavemin=0., wavemax=60000., Smin=0.):
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

    for line in filehandle:

        if (line[0] == '#'):
            continue

        if (len(line) < 160):
            raise ImportError('The imported file ("' + filename + '") does not appear to be a HITRAN2012-format data file.')

        vc = np.float64(line[3:15])
        S = np.float64(line[15:25])
        if vc > wavemax:
            break  # don't need to continue to save
        if ((wavemin <= vc) and (S > Smin)):
            linelist['M'].append(np.uint(line[0:2]))
            I = line[2]
            if I == '0':
                linelist['I'].append(np.uint(10))
            elif I == 'a':
                linelist['I'].append(np.uint(11))
            elif I == 'b':
                linelist['I'].append(np.uint(11))
            else:
                linelist['I'].append(np.uint(line[2]))
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
            linelist['Ierr'].append(np.uint(line[127:133]))
            linelist['Iref'].append(np.uint(line[133:145]))
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


def calculate_hitran_xsec(linelist, M, wavenumbers, T=296.e0, P=101325., Pref=101325., qmix=0.0, wreach=25.0):
    """
    Given the HITRAN linelist (line centers and line strengths) for a molecule, digitize the result into a spectrum of
    absorption cross-section in units of cm^2.

    Parameters
    ----------
    linelist : dict of ndarrays
        The HITRAN data corresponding to a given molecule.
    M : molecule number
    wavemin : float, optional
        The minimum wavenumber os the spectral region of interest.
    wavemax : float, optional
        The maximum wavenumber os the spectral region of interest.
    units : str, optional
        A string describing in what units of the output cross-section should be given in. Choices available are:
        {'cm^2/mole', 'cm^2.ppm', 'm^2/mole', 'm^2.ppm', 'm^2', cm^2}.

    T : float
        Temperature in Kelvin (Tref = 296. K)
    P : float
        Pressure in Pascal (Pref = 1 atm = 101325 Pa)

    Returns
    -------
    xsec : array_like
        The mean absorption cross-section (in cm^2) per molecule, evaluated at the wavenumbers given by input `waves`.
    """

    import scipy.constants as constants

    kb = constants.value('Boltzmann constant')*1e7
    c = constants.value('speed of light in vacuum')*1e2
    amu = constants.value('atomic mass constant')*1e3
    h = constants.value('Planck constant')*1e7
    hck = h*c/kb

    #print(kb, c, amu, h, hck)
   
    xsec = np.zeros_like(wavenumbers)

    # check for errors in Epp (was an issue in previous versions of HITRAN)
    if any(linelist['Epp']<0.0):
        raise ValueError("Need to check for error flags in Epp")

    # scale line strengths to temperature

    #print(hck, np.max(linelist['S']))
    hckt = hck*(1./T-1./296.)
    qratio = [qtips(296.0,M,i)/qtips(T,M,i) for i in range(1,get_molecule_nisops(M)+1)]
    #linestrengths = linestrengths*(qtips(296., M)/qtips(T, M)) *exp(-hckt*Epps)
    #linestrengths = [S*qratio[int(I-1)]*exp(-hckt*Epp) for (S,I,Epp) in zip(linestrengths,isops,Epps)]
    linestrengths = [S*qratio[int(I-1)]*np.exp(-hckt*Epp)*(1.0-np.exp(-hck*vc/T))/(1.0-np.exp(-hck*vc/296.0))
                     for (vc,S,I,Epp) in zip(linelist['vc'],linelist['S'],linelist['I'],linelist['Epp'])]

    # pre-calculate doppler widths and lorentz widths
    vcs = linelist['vc'] + linelist['delta']*(P/Pref)
    masses = [get_iso_mass(M,i)*amu for i in range(1,get_molecule_nisops(M)+1)]
    alphas = [np.sqrt(2.0*kb*T/masses[int(I-1)]) for I in linelist['I']] * vcs/c
    gammas = ((1.0-qmix)*linelist['gamma-air']+qmix*linelist['gamma-self']) * (P/Pref) * (296./T)**linelist['N']

    #print(P, Pref, T, 296.)
    #for il, linecenter, linestrength, alpha, gamma in zip(range(len(vcs)), vcs, linestrengths, alphas, gammas):
        #print(linecenter, linestrength, alpha/np.sqrt(2.)*np.sqrt(2.*np.log(2.)), gamma) 
        #print(linelist['vc'][il], linelist['S'][il], linelist['gamma-air'][il], linelist['gamma-self'][il]) 
        #print((1.0-qmix)*linelist['gamma-air'][il], qmix*linelist['gamma-self'][il], (296./T)**linelist['N'][il])
        #print(qratio[0], masses[0], np.exp(-hckt*linelist['Epp'][il]), (1.0-np.exp(-hck*linelist['vc'][il]/T))/(1.0-np.exp(-hck*linelist['vc'][il]/296.0)))

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



### ======================================================
#def downsample_spectrum(waves, spectrum, downsampled_waves=None, downsampled_channel_boundaries=None):
#    '''
#    (Right now, we assume uniformly sampled wavelengths/wavenumbers.)
#    '''
#
#    nwaves = np.alen(waves)
#
#    ## If it is not already defined, make the list of channel boundary wavelengths.
#    if (downsampled_waves != None) and (downsampled_channel_boundaries == None):
#        dw = downsampled_waves[1] - downsampled_waves[0]
#        downsampled_channel_boundaries = append(amin(downsampled_waves)-(dw/2.0), downsampled_waves+(dw/2.0))
#    elif (downsampled_waves == None) and (downsampled_channel_boundaries == None):
#        raise ValueError('Either "downsampled_waves" or "downsampled_channel_boundaries" is required as an input.')
#
#    ## Generate the channel basis functions used to represent the low-resolution spectral channels in terms
#    ## of the high-resolution data.
#    nchannels = alen(downsampled_channel_boundaries) - 1
#    #print('downwaves=', downwaves)
#    #print('downsampled_channel_boundaries=', downsampled_channel_boundaries)
#    downspectrum = zeros(nchannels)
#
#    ## From the list of channel boundary wavelengths, generate the channel basis functions.
#    for n in arange(nchannels):
#        okay = (waves > downsampled_channel_boundaries[n]) * (waves <= downsampled_channel_boundaries[n+1])
#        downspectrum[n] = nanmean(spectrum[okay])
#
#    return(downsampled_channel_boundaries, downspectrum)
#
#
### =================================================================================================
#def draw_block_spectrum(channel_boundaries, spectrum, newfigure=True, title=None, **kwargs):
#    '''
#    Draw a plot where the spectral channels are nonuniform in width and shaped by histogram-like rectangles.
#
#    Parameters
#    ----------
#    channel_boundaries : array_like
#        A vector of length Nw+1 giving the wavelengths defining the boundaries of the Nw spectral channels.
#    spectrum : array_like
#        A Nw length vector.
#    newfigure : bool
#        Whether or not to call matplotlib's `figure()` function.
#    title : str
#        The plot title.
#    **kwargs : any
#        Any keyword arguments that can be used by matplotlib's `plot()` function.
#    '''
#
#    from matplotlib import pyplot as plt
#
#    channel_boundaries = asarray(channel_boundaries)
#    spectrum = asarray(spectrum)
#    assert (np.alen(channel_boundaries) == 1 + np.alen(spectrum)), 'Input "channel_boundaries" must have length 1 more than input "spectrum".'
#
#    cb = channel_boundaries
#    nchannels = alen(cb) - 1
#
#    x = []
#    y = []
#    x.append(cb[0])
#    y.append(0.0)
#
#    for n in arange(nchannels):
#        x.append(cb[n])
#        x.append(cb[n+1])
#        y.append(spectrum[n])
#        y.append(spectrum[n])
#
#    x.append(cb[-1])
#    y.append(0.0)
#
#    if newfigure:
#        fig = plt.figure()
#        if (title != None): fig.canvas.set_window_title(title)
#        xmin = amin(x)
#        xmax = amax(x)
#        xptp = xmax - xmin
#        xmean = 0.5 * (xmin + xmax)
#        xlo = xmean - 0.55 * xptp
#        xhi = xmean + 0.55 * xptp
#
#        ymin = amin(y)
#        ymax = amax(y)
#        yptp = ymax - ymin
#        ymean = 0.5 * (ymin + ymax)
#        ylo = ymean - 0.55 * yptp
#        yhi = ymean + 0.55 * yptp
#    else:
#        ## First grab the existing plot limits. If these were previously determined by draw_block_spectrum(),
#        ## then we need only rescale the plot range by (0.5/0.55) to get to the original data limits.
#        ## Appending these to the current data vector, we can update the plot limits using both old and new
#        ## data. First grab the existing limits.
#        (x0,x1,y0,y1) = plt.axis()
#
#        x0mean = 0.5 * (x0 + x1)
#        x0ptp = (0.5 / 0.55) * (x1 - x0)
#        x0min = x0mean - 0.5 * x0ptp
#        x0max = x0mean + 0.5 * x0ptp
#
#        y0mean = 0.5 * (y0 + y1)
#        y0ptp = (0.5 / 0.55) * (y1 - y0)
#        y0min = y0mean - 0.5 * y0ptp
#        y0max = y0mean + 0.5 * y0ptp
#
#        ## Next determine the new plot range using the old and new data limits.
#        xmin = amin(append(array(x), x0min))
#        xmax = amax(append(array(x), x0max))
#        xptp = xmax - xmin
#        xmean = 0.5 * (xmin + xmax)
#        xlo = xmean - 0.55 * xptp
#        xhi = xmean + 0.55 * xptp
#
#        ymin = amin(append(array(y), y0min))
#        ymax = amax(append(array(y), y0max))
#        yptp = ymax - ymin
#        ymean = 0.5 * (ymin + ymax)
#        ylo = ymean - 0.55 * yptp
#        yhi = ymean + 0.55 * yptp
#
#    plt.plot(x, y, **kwargs)
#    if (title != None): plt.title(title)
#    if newfigure:
#        plt.axis([xlo,xhi,ylo,yhi])
#
#    return


