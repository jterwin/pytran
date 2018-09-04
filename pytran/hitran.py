

from __future__ import print_function, division
from numpy import *
import matplotlib.pyplot as plt
import matplotlib.cm as cm

from pytran.hitran_utils import *


__all__ = ['lorentzian_profile', 'doppler_profile', 'voigt_profile', 'read_hitran2012_parfile',
               'calculate_hitran_xsec']

## ======================================================
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
    L = (S / pi) * gamma / ((wavn - wavn0)**2 + gamma**2)
    return(L)
    
## =========
def doppler_profile(wavn, S, alpha, wavn0):
    '''
    doppler profile
    '''
    G = (S / (sqrt(pi)*alpha)) * exp(-((wavn-wavn0) / alpha)**2)
    return(G)
    
## =========
def voigt_profile(wavn, S, alpha, gamma, wavn0):
    '''
    voigt profile
    '''
    from scipy.special import wofz
    z = ((wavn-wavn0) + 1j*gamma)/(alpha)
    V = S/(sqrt(pi)*alpha)*wofz(z).real
    return(V)

## ======================================================
def read_hitran2012_parfile(filename, wavemin=0., wavemax=60000.):
    """
    Given a HITRAN2012-format text file, read in the parameters of the molecular absorption features.

    Parameters
    ----------
    filename : str
        The filename to read in.

    Return
    ------
    data : dict
        The dictionary of HITRAN data for the molecule.
    """

    if filename.endswith('.zip'):
        import zipfile
        zip = zipfile.ZipFile(filename, 'r')
        (object_name, ext) = os.path.splitext(os.path.basename(filename))
        print(object_name, ext)
        filehandle = zip.read(object_name).splitlines()
    else:
        filehandle = open(filename, 'r')

    data = {'M':[],               ## molecule identification number
            'I':[],               ## isotope number
            'linecenter':[],      ## line center wavenumber (in cm^{-1})
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
            'gpp':[]}             ## statistical weight of the lower state

    print('Reading "' + filename + '" ...')

    for line in filehandle:
        if (len(line) < 160):
            raise ImportError, 'The imported file ("' + filename + '") does not appear to be a HITRAN2012-format data file.'

        if (wavemin <= float64(line[3:15]) <= wavemax):
            data['M'].append(uint(line[0:2]))
            I = line[2]
            if I == '0':
                data['I'].append(uint(10))
            elif I == 'a':
                data['I'].append(uint(11))
            elif I == 'b':
                data['I'].append(uint(11))
            else:
                data['I'].append(uint(line[2]))
            data['linecenter'].append(float64(line[3:15]))
            data['S'].append(float64(line[15:25]))
            data['Acoeff'].append(float64(line[25:35]))
            data['gamma-air'].append(float64(line[35:40]))
            data['gamma-self'].append(float64(line[40:45]))
            data['Epp'].append(float64(line[45:55]))
            data['N'].append(float64(line[55:59]))
            data['delta'].append(float64(line[59:67]))
            data['Vp'].append(line[67:82])
            data['Vpp'].append(line[82:97])
            data['Qp'].append(line[97:112])
            data['Qpp'].append(line[112:127])
            data['Ierr'].append(line[127:133])
            data['Iref'].append(line[133:145])
            data['flag'].append(line[145])
            data['gp'].append(line[146:153])
            data['gpp'].append(line[153:160])

    if filename.endswith('.zip'):
        zip.close()
    else:
        filehandle.close()

    for key in data:
        data[key] = array(data[key])

    return data



    
## ======================================================
def calculate_hitran_xsec(data, M, wnmin=None, wnmax=None, npts=20001, T=296.e0, p=1.0e6, units='cm^2'):
    """
    Given the HITRAN data (line centers and line strengths) for a molecule, digitize the result into a spectrum of
    absorption cross-section in units of cm^2.

    Parameters
    ----------
    data : dict of ndarrays
        The HITRAN data corresponding to a given molecule.
    M : molecule number
    wavemin : float, optional
        The minimum wavenumber os the spectral region of interest.
    wavemax : float, optional
        The maximum wavenumber os the spectral region of interest.
    units : str, optional
        A string describing in what units of the output cross-section should be given in. Choices available are:
        {'cm^2/mole', 'cm^2.ppm', 'm^2/mole', 'm^2.ppm', 'm^2', cm^2}.

    Returns
    -------
    waves : ndarray
        The wavenumbers at which the cross-section data is evaluated.
    xsec : array_like
        The mean absorption cross-section (in cm^2) per molecule, evaluated at the wavenumbers given by input `waves`.
    """
    
    import scipy.constants as constants
    kb = constants.value('Boltzmann constant')*1e7
    c = constants.value('speed of light in vacuum')*1e2
    amu = constants.value('atomic mass constant')*1e3
    h = constants.value('Planck constant')*1e7
    
    if (wnmin == None):
        wnmin = amin( data['linecenter']) - 0.1
    if (wnmax == None):
        wnmax = amax( data['linecenter']) + 0.1

    ## First step: 
    okay = (data['I'] >= 1)
    #okay = ((wavemin-10.) <= data['linecenter'] <= (wavemax+10.))
    #okay = (data['I'] == 1) * (data['linecenter']>=(wavemin-10.)) * (data['linecenter']<=wavemax+10.)

    # pull some arrays out of line list
    linecenters = array(data['linecenter'][okay])       ## line centers in wavenumbers
    linestrengths = array(data['S'][okay])
    Epps = array(data['Epp'][okay])
    isops = array(data['I'][okay])
    gammas = array(data['gamma-air'][okay])
    Ns = array(data['N'])[okay]
    
    # check for errors in Epp (was an issue in previous versions of HITRAN)
    if any(Epps<0.0):
        print("need to check for error flags in Epp")
        import sys
        sys.exit()

    # scale line strengths to temperature
    hck = h*c/kb
    hckt = hck*(1./T-1./296.)
    qratio = [qtips(296.0,M,i)/qtips(T,M,i) for i in range(1,get_molecule_nisops(M)+1)]
    #linestrengths = linestrengths*(qtips(296., M)/qtips(T, M)) *exp(-hckt*Epps)
    #linestrengths = [S*qratio[int(I-1)]*exp(-hckt*Epp) for (S,I,Epp) in zip(linestrengths,isops,Epps)]
    linestrengths = [S*qratio[int(I-1)]*exp(-hckt*Epp)*(1.0-exp(-hck*vc/T))/(1.0-exp(-hck*vc/296.0))
                     for (vc,S,I,Epp) in zip(linecenters,linestrengths,isops,Epps)]

    # pre-calculate doppler widths and lorentz widths
    #mass = get_molecule_mass(M)*amu
    #alphas = ones_like(linecenters)*sqrt(2.0*kb*T/mass)
    masses = [get_iso_mass(M,i) for i in range(1,get_molecule_nisops(M)+1)]
    alphas = [sqrt(2.0*kb*T/masses[int(I-1)]) for I in isops]
    gammas = gammas*(p/1e6)/(T/296.)**Ns
    
    ## Create a spectrum linearly sampled in wavenumber
    wavenumbers = linspace(wnmin, wnmax, npts)
    xsec = zeros_like(wavenumbers)

    nlines = alen(linecenters)
    for i in arange(nlines):
        linecenter = linecenters[i]
        linestrength = linestrengths[i]
        gamma = gammas[i]
        alpha = alphas[i]*linecenter/c

        ## If the spectral line is well outside our region of interest, then ignore it.
        #if (linecenter < amin(wavenumbers-0.5)):
        #    continue
        #elif (linecenter > amax(wavenumbers+0.5)):
        #    continue

        ## Note: the quantity sum(L * dk) should sum to "S"!
        #L = lorentzian_profile(wavenumbers, linestrength, gamma, linecenter)
        #xsec += L
        #D = doppler_profile(wavenumbers, linestrength, alpha, linecenter)
        #xsec += D
        V = voigt_profile(wavenumbers, linestrength, alpha, gamma, linecenter)
        xsec += V

    if units.endswith('/mole'):
        xsec = xsec * 6.022E23
    elif units.endswith('.ppm'):
        xsec = xsec * 2.686E19

    if units.startswith('cm^2'):
        pass
    elif units.startswith('m^2'):
        xsec = xsec / 10000.0

    return(wavenumbers, xsec)



## ======================================================
def downsample_spectrum(waves, spectrum, downsampled_waves=None, downsampled_channel_boundaries=None):
    '''
    (Right now, we assume uniformly sampled wavelengths/wavenumbers.)
    '''

    nwaves = alen(waves)

    ## If it is not already defined, make the list of channel boundary wavelengths.
    if (downsampled_waves != None) and (downsampled_channel_boundaries == None):
        dw = downsampled_waves[1] - downsampled_waves[0]
        downsampled_channel_boundaries = append(amin(downsampled_waves)-(dw/2.0), downsampled_waves+(dw/2.0))
    elif (downsampled_waves == None) and (downsampled_channel_boundaries == None):
        raise ValueError, 'Either "downsampled_waves" or "downsampled_channel_boundaries" is required as an input.'

    ## Generate the channel basis functions used to represent the low-resolution spectral channels in terms
    ## of the high-resolution data.
    nchannels = alen(downsampled_channel_boundaries) - 1
    #print('downwaves=', downwaves)
    #print('downsampled_channel_boundaries=', downsampled_channel_boundaries)
    downspectrum = zeros(nchannels)

    ## From the list of channel boundary wavelengths, generate the channel basis functions.
    for n in arange(nchannels):
        okay = (waves > downsampled_channel_boundaries[n]) * (waves <= downsampled_channel_boundaries[n+1])
        downspectrum[n] = nanmean(spectrum[okay])

    return(downsampled_channel_boundaries, downspectrum)

## =================================================================================================
def draw_block_spectrum(channel_boundaries, spectrum, newfigure=True, title=None, **kwargs):
    '''
    Draw a plot where the spectral channels are nonuniform in width and shaped by histogram-like rectangles.

    Parameters
    ----------
    channel_boundaries : array_like
        A vector of length Nw+1 giving the wavelengths defining the boundaries of the Nw spectral channels.
    spectrum : array_like
        A Nw length vector.
    newfigure : bool
        Whether or not to call matplotlib's `figure()` function.
    title : str
        The plot title.
    **kwargs : any
        Any keyword arguments that can be used by matplotlib's `plot()` function.
    '''

    channel_boundaries = asarray(channel_boundaries)
    spectrum = asarray(spectrum)
    assert (alen(channel_boundaries) == 1 + alen(spectrum)), 'Input "channel_boundaries" must have length 1 more than input "spectrum".'

    cb = channel_boundaries
    nchannels = alen(cb) - 1

    x = []
    y = []
    x.append(cb[0])
    y.append(0.0)

    for n in arange(nchannels):
        x.append(cb[n])
        x.append(cb[n+1])
        y.append(spectrum[n])
        y.append(spectrum[n])

    x.append(cb[-1])
    y.append(0.0)

    if newfigure:
        fig = plt.figure()
        if (title != None): fig.canvas.set_window_title(title)
        xmin = amin(x)
        xmax = amax(x)
        xptp = xmax - xmin
        xmean = 0.5 * (xmin + xmax)
        xlo = xmean - 0.55 * xptp
        xhi = xmean + 0.55 * xptp

        ymin = amin(y)
        ymax = amax(y)
        yptp = ymax - ymin
        ymean = 0.5 * (ymin + ymax)
        ylo = ymean - 0.55 * yptp
        yhi = ymean + 0.55 * yptp
    else:
        ## First grab the existing plot limits. If these were previously determined by draw_block_spectrum(),
        ## then we need only rescale the plot range by (0.5/0.55) to get to the original data limits.
        ## Appending these to the current data vector, we can update the plot limits using both old and new
        ## data. First grab the existing limits.
        (x0,x1,y0,y1) = plt.axis()

        x0mean = 0.5 * (x0 + x1)
        x0ptp = (0.5 / 0.55) * (x1 - x0)
        x0min = x0mean - 0.5 * x0ptp
        x0max = x0mean + 0.5 * x0ptp

        y0mean = 0.5 * (y0 + y1)
        y0ptp = (0.5 / 0.55) * (y1 - y0)
        y0min = y0mean - 0.5 * y0ptp
        y0max = y0mean + 0.5 * y0ptp

        ## Next determine the new plot range using the old and new data limits.
        xmin = amin(append(array(x), x0min))
        xmax = amax(append(array(x), x0max))
        xptp = xmax - xmin
        xmean = 0.5 * (xmin + xmax)
        xlo = xmean - 0.55 * xptp
        xhi = xmean + 0.55 * xptp

        ymin = amin(append(array(y), y0min))
        ymax = amax(append(array(y), y0max))
        yptp = ymax - ymin
        ymean = 0.5 * (ymin + ymax)
        ylo = ymean - 0.55 * yptp
        yhi = ymean + 0.55 * yptp

    plt.plot(x, y, **kwargs)
    if (title != None): plt.title(title)
    if newfigure:
        plt.axis([xlo,xhi,ylo,yhi])

    return


