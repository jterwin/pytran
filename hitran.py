

from __future__ import print_function, division
from numpy import *
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import os
import glob

__all__ = ['lorentzian_profile', 'read_hitran2012_parfile', 'translate_molecule_identifier',
           'get_molecule_identifier', 'calculate_hitran_xsec', 'downsample_spectrum', 'draw_block_spectrum']

## ======================================================
def lorentzian_profile(kappa, S, gamma, kappa0):
    '''
    Calculate a Lorentzian absorption profile.

    Parameters
    ----------
    kappa : ndarray
        The array of wavenumbers at which to sample the profile.
    S : float
        The absorption line "strength" (the integral of the entire line is equal to S).
    gamma : float
        The linewidth parameter of the profile.
    kappa0 : float
        The center position of the profile (in wavenumbers).

    Returns
    -------
    L : ndarray
        The sampled absorption profile.
    '''
    L = (S / pi) * gamma / ((kappa - kappa0)**2 + gamma**2)
    return(L)
    
## =========
def doppler_profile(wavn, S, alpha, wav0):
    '''
    doppler profile
    '''
    G = (S / (sqrt(pi)*alpha)) * exp(-((wavn-wav0) / alpha)**2)
    return(G)
    
## =========
def voigt_profile(wavn, S, alpha, gamma, wav0):
    '''
    voigt profile
    '''
    from scipy.special import wofz
    z = ((wavn-wav0) + 1j*gamma)/(alpha)
    V = S/(sqrt(pi)*alpha)*wofz(z).real
    return(V)

## ======================================================
def read_hitran2012_parfile(filename, wavemin=0., wavemax=60000.):
    '''
    Given a HITRAN2012-format text file, read in the parameters of the molecular absorption features.

    Parameters
    ----------
    filename : str
        The filename to read in.

    Return
    ------
    data : dict
        The dictionary of HITRAN data for the molecule.
    '''

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

    return(data)

## ======================================================
def translate_molecule_identifier(M):
    '''
    For a given input molecule identifier number, return the corresponding molecular formula.

    Parameters
    ----------
    M : int
        The HITRAN molecule identifier number.

    Returns
    -------
    molecular_formula : str
        The string describing the molecule.
    '''

    trans = { '1':'H2O',    '2':'CO2',   '3':'O3',      '4':'N2O',   '5':'CO',    '6':'CH4',   '7':'O2',     '8':'NO',
              '9':'SO2',   '10':'NO2',  '11':'NH3',    '12':'HNO3', '13':'OH',   '14':'HF',   '15':'HCl',   '16':'HBr',
             '17':'HI',    '18':'ClO',  '19':'OCS',    '20':'H2CO', '21':'HOCl', '22':'N2',   '23':'HCN',   '24':'CH3Cl',
             '25':'H2O2',  '26':'C2H2', '27':'C2H6',   '28':'PH3',  '29':'COF2', '30':'SF6',  '31':'H2S',   '32':'HCOOH',
             '33':'HO2',   '34':'O',    '35':'ClONO2', '36':'NO+',  '37':'HOBr', '38':'C2H4', '39':'CH3OH', '40':'CH3Br',
             '41':'CH3CN', '42':'CF4',  '43':'C4H2',   '44':'HC3N', '45':'H2',   '46':'CS',   '47':'SO3'}
    return(trans[str(M)])

## ======================================================
def get_molecule_identifier(molecule_name):
    '''
    For a given input molecular formula, return the corresponding HITRAN molecule identifier number.

    Parameters
    ----------
    molecular_formula : str
        The string describing the molecule.

    Returns
    -------
    M : int
        The HITRAN molecular identified number.
    '''

    trans = { '1':'H2O',    '2':'CO2',   '3':'O3',      '4':'N2O',   '5':'CO',    '6':'CH4',   '7':'O2',     '8':'NO',
              '9':'SO2',   '10':'NO2',  '11':'NH3',    '12':'HNO3', '13':'OH',   '14':'HF',   '15':'HCl',   '16':'HBr',
             '17':'HI',    '18':'ClO',  '19':'OCS',    '20':'H2CO', '21':'HOCl', '22':'N2',   '23':'HCN',   '24':'CH3Cl',
             '25':'H2O2',  '26':'C2H2', '27':'C2H6',   '28':'PH3',  '29':'COF2', '30':'SF6',  '31':'H2S',   '32':'HCOOH',
             '33':'HO2',   '34':'O',    '35':'ClONO2', '36':'NO+',  '37':'HOBr', '38':'C2H4', '39':'CH3OH', '40':'CH3Br',
             '41':'CH3CN', '42':'CF4',  '43':'C4H2',   '44':'HC3N', '45':'H2',   '46':'CS',   '47':'SO3'}
    ## Invert the dictionary.
    trans = {v:k for k,v in trans.items()}
    return(int(trans[molecule_name]))
    
## ======================================================
def molecule_mass(M):
    
    masses = { '1':18.01,    '2':43.99,   '3':47.98,     '4':44.00,   '5':27.99,   '6':16.03,   '7':31.99,    '8':30.00,
               '9':63.96,   '10':45.99,  '11':17.03,    '12':63.00,  '13':17.00,  '14':20.01,  '15':35.98,   '16':79.93,
              '17':127.91,  '18':50.96,  '19':59.97,    '20':30.01,  '21':51.97,  '22':28.01,  '23':27.01,   '24':50.00,
              '25':34.01,   '26':26.02,  '27':30.05,    '28':34.00,  '29':65.99,  '30':145.96, '31':33.99,   '32':46.01,
              '33':33.00,   '34':15.99,  '35':96.96,    '36':30.00,  '37':95.92,  '38':28.03,  '39':32.03,   '40':93.94,
              '41':41.03,   '42':87.99,  '43':50.00,    '44':61.00,  '45':52.00,  '46':43.97,  '47':2.01,    '48':48.00,
              '49':40.00,   '50':15.00,  '51':76.00}
    return(masses[str(M)])



Tdata = [60.,  85., 110., 135., 160., 185., 210., 235.,
         260., 285., 310., 335., 360., 385., 410., 435., 460., 485.,
         510., 535., 560., 585., 610., 635., 660., 685., 710., 735.,
         760., 785., 810., 835., 860., 885., 910., 935., 960., 985.,
         1010.,1035.,1060.,1085.,1110.,1135.,1160.,1185.,1210.,1235.,
         1260.,1285.,1310.,1335.,1360.,1385.,1410.,1435.,1460.,1485.,
         1510.,1535.,1560.,1585.,1610.,1635.,1660.,1685.,1710.,1735.,
         1760.,1785.,1810.,1835.,1860.,1885.,1910.,1935.,1960.,1985.,
         2010.,2035.,2060.,2085.,2110.,2135.,2160.,2185.,2210.,2235.,
         2260.,2285.,2310.,2335.,2360.,2385.,2410.,2435.,2460.,2485.,
         2510.,2535.,2560.,2585.,2610.,2635.,2660.,2685.,2710.,2735.,
         2760.,2785.,2810.,2835.,2860.,2885.,2910.,2935.,2960.,2985.,
         3010.]

qt_tab = [0.54800E+02, 0.91500E+02, 0.13410E+03,
          0.18180E+03, 0.23410E+03, 0.29070E+03, 0.35140E+03, 0.41600E+03,
          0.48450E+03, 0.55720E+03, 0.63420E+03, 0.71600E+03, 0.80310E+03,
          0.89590E+03, 0.99520E+03, 0.11017E+04, 0.12161E+04, 0.13393E+04,
          0.14721E+04, 0.16155E+04, 0.17706E+04, 0.19384E+04, 0.21202E+04,
          0.23172E+04, 0.25307E+04, 0.27624E+04, 0.30137E+04, 0.32864E+04,
          0.35823E+04, 0.39034E+04, 0.42519E+04, 0.46300E+04, 0.50402E+04,
          0.54853E+04, 0.59679E+04, 0.64913E+04, 0.70588E+04, 0.76739E+04,
          0.83404E+04, 0.90625E+04, 0.98446E+04, 0.10691E+05, 0.11608E+05,
          0.12600E+05, 0.13674E+05, 0.14835E+05, 0.16090E+05, 0.17447E+05,
          0.18914E+05, 0.20500E+05, 0.22212E+05, 0.24063E+05, 0.26061E+05,
          0.28218E+05, 0.30548E+05, 0.33063E+05, 0.35778E+05, 0.38708E+05,
          0.41871E+05, 0.45284E+05, 0.48970E+05, 0.52940E+05, 0.57230E+05,
          0.61860E+05, 0.66860E+05, 0.72250E+05, 0.78070E+05, 0.84350E+05,
          0.91130E+05, 0.98450E+05, 0.10635E+06, 0.11488E+06, 0.12408E+06,
          0.13403E+06, 0.14480E+06, 0.15640E+06, 0.16890E+06, 0.18240E+06,
          0.19700E+06, 0.21280E+06, 0.22980E+06, 0.24830E+06, 0.26820E+06,
          0.28970E+06, 0.31290E+06, 0.33800E+06, 0.36520E+06, 0.39450E+06,
          0.42600E+06, 0.46000E+06, 0.49700E+06, 0.53700E+06, 0.58100E+06,
          0.62700E+06, 0.67800E+06, 0.73300E+06, 0.79200E+06, 0.85600E+06,
          0.92500E+06, 0.10000E+07, 0.10800E+07, 0.11670E+07, 0.12610E+07,
          0.13620E+07, 0.14720E+07, 0.15910E+07, 0.17190E+07, 0.18600E+07,
          0.20100E+07, 0.21700E+07, 0.23400E+07, 0.25300E+07, 0.27300E+07,
          0.29500E+07, 0.31800E+07, 0.34300E+07, 0.37000E+07, 0.39900E+07,
          0.42856E+07]




    
## ======================================================
def calculate_hitran_xsec(data, M, wnmin=None, wnmax=None, npts=20001, T=296.e0, p=1.0e6, units='m^2'):
    '''
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
    '''
    
    import scipy.constants as constants
    kb = constants.physical_constants['Boltzmann constant'][0]*1e7
    c = constants.physical_constants['speed of light in vacuum'][0]*1e2
    amu = constants.physical_constants['atomic mass constant'][0]*1e3
    h = constants.physical_constants['Planck constant'][0]*1e7
    

    # TIPS function for methane
    Z = polyfit(Tdata,log(qt_tab),10)
    ZZ = poly1d(Z)
    qt_fun = lambda x: exp(ZZ(x))
    print("%f %f %f %f %f" % (qt_fun(40.),qt_fun(60.),qt_fun(80.),qt_fun(100.),qt_fun(0.)) )

    
    if (wnmin == None):
        wnmin = amin( data['linecenter']) - 0.1
    if (wnmax == None):
        wnmax = amax( data['linecenter']) + 0.1

    ## First step: remove any data points that do not correspond to the primary isotope. (If we want to use isotopes,
    ## then we need to figure out mixing ratios.) For most HITRAN gases, the primary isotope is about 99% of the total
    ## atmospheric composition.
    okay = (data['I'] >= 1)
    #okay = ((wavemin-10.) <= data['linecenter'] <= (wavemax+10.))
    #print(data['I'])
    #print(okay)
    #okay = (data['I'] == 1) * (data['linecenter']>=(wavemin-10.)) * (data['linecenter']<=wavemax+10.)
    #print("%d,%d" % (len(data['S']),len(okay)) )

    linecenters = array(data['linecenter'][okay])       ## line centers in wavenumbers
    linestrengths = array(data['S'][okay])
    Epp = array(data['Epp'][okay])

    if any(Epp<0.0):
        print("need to check for error flags in Epp")
        import sys
        sys.exit()

    hckt = h*c/kb*(1./T-1./296.)
    linestrengths = linestrengths*(qt_fun(296.)/qt_fun(T)) *exp(-hckt*Epp)
    #print("%e, %e, %e, %e" % (h*c/kb,hckt,amin(Epp),amax(Epp)) )

    
    mass = molecule_mass(M)*amu
    alphas = ones_like(linecenters)*sqrt(2.0*kb*T/mass)
    
    gammas = array(data['gamma-air'][okay])
    Ns = array(data['N'])[okay]
    gammas = gammas*(p/1e6)/(T/296.)**Ns
    
    #print(mass,alphas[0]*linecenters[0]/c,gammas[0])

    

    ## Create a spectrum linearly sampled in wavenumber
    wavenumbers = linspace(wnmin, wnmax, npts)
    # waves = 10000.0 / wavenumbers
    xsec = zeros_like(wavenumbers)

    ## Define the list of channel boundary wavelengths.
    #dk = wavenumbers[1] - wavenumbers[0]

    nlines = alen(linecenters)
    print("%d"%nlines)

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
        
        #print(i, mass, alpha, gamma , linecenter, alphas[i])
        #import sys
        #sys.exit()
        
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


## ==================================================================================================
## ==================================================================================================

if (__name__ == "__main__"):
    #molecule = 'H2'       ## water
    #molecule = 'CO2'       ## carbon dioxide
    #molecule = 'NH3'       ## ammonia
    #molecule = 'SO2'       ## sulfur dioxide
    molecule = 'CH4'       ## methane
    #molecule = 'H2S'       ## hydrogen sulfide
    #molecule = 'O3'        ## ozone
    #molecule = 'C2H6'      ## ethane

    #units = 'm^2/mole'
    #units = 'm^2.ppm'
    #units = 'cm^2/mole'
    #units = 'cm^2.ppm'
    #units = 'm^2'
    units = 'cm^2'

    #wnmin, wnmax = (4000., 4600.)
    #wnmin, wnmax = (2800., 3200.)
    wnmin, wnmax = (1090., 1225.)
    #wnmin, wnmax = (1225., 1360.)
    print('(%f.1, %f.1)' % (wnmin,wnmax))

    
    M = get_molecule_identifier(molecule)
    filenames = glob.glob('../hitran/HITRAN2012/By-Molecule/Uncompressed-files/%02i_hit12.*' % M)
    filename = filenames[0]

    filename = os.path.normpath(os.path.abspath(filename))
    data = read_hitran2012_parfile(filename, wnmin, wnmax)
    nlines = len(data['S'])
    print('Found %i lines' % nlines)

    
    ## Next set up for calculating the absorption cross-section, given the HITRAN database's values for the
    ## line centers and line strengths.
    print('Calculating the absorption cross-section spectrum ...')
    (wns, xsec) = calculate_hitran_xsec(data, M, wnmin, wnmax, units=units, T=100.,p=10e1)

     
    print("plotting")
    fig, ax = plt.subplots()
    fig.canvas.set_window_title(molecule)
    ax.semilogy(wns, xsec, 'k-')
    
    ax.set_title(molecule)
    ax.set_ylabel('Cross-section (' + units + ')')
    ax.set_xlabel('wavenumber (cm$^{-1}$)')
    ax.set_xlim([wnmin-0.05*(wnmax-wnmin),wnmax+0.05*(wnmax-wnmin)])



    
    print("sorting")
    fig, ax = plt.subplots()
    ind = xsec.argsort()
    ax.semilogy(linspace(0,1,xsec.size),xsec[ind])


    sigmabar = sum(xsec)
    print( '%e' % (sigmabar) )
    sigmabar = sigmabar*(wns[1]-wns[0])  
    print ( '%e' % (sigmabar) )
    print ( '%e' % (sigmabar/(wns[-1]-wns[0])) )
    plt.show()
