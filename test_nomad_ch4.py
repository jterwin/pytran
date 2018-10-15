


## ==================================================================================================
## ==================================================================================================

if (__name__ == "__main__"):

    import os
    import glob
    
    import pytran
    import numpy as np
    from matplotlib import pyplot as plt

    import h5py

    #molecule = 'H2'       ## water
    #molecule = 'CO2'       ## carbon dioxide
    #molecule = 'NH3'       ## ammonia
    #molecule = 'SO2'       ## sulfur dioxide
    molecule = 'CH4'       ## methane
    #molecule = 'H2S'       ## hydrogen sulfide
    #molecule = 'O3'        ## ozone
    #molecule = 'C2H6'      ## ethane

    units = 'cm^2'

    wnmin, wnmax = (3005., 3040.)
    print('(%.1f, %.1f)' % (wnmin,wnmax))

    #HITRANDIR = os.path.expanduser('~/lib/HITRAN/HITRAN2012/By-Molecule/Uncompressed-files')
    HITRANDIR = os.path.expanduser('~/projects/NOMAD/Science/Radiative_Transfer/Auxiliary_files/Spectroscopy')
    print(HITRANDIR)

    M = pytran.get_molecule_id(molecule)    
    filenames = glob.glob(os.path.join(HITRANDIR,'%02i_hit*.par' % M))
    filename = filenames[0]
    filename = os.path.normpath(os.path.abspath(filename))
    print(filename)
    
    data = pytran.read_hitran2012_parfile(filename, wnmin, wnmax, Smin=1.e-25)
    nlines = len(data['S'])
    print('Found %i lines' % nlines)


    ## Next set up for calculating the absorption cross-section, given the HITRAN database's values for the
    ## line centers and line strengths.
    print('Calculating the absorption cross-section spectrum ...')
    dwn = 0.001
    Nbwn = np.ceil((wnmax-wnmin)/dwn) + 1
    print(Nbwn)
    wns = np.linspace(wnmin, wnmax, Nbwn)

    if False:    # plot single spectrum
        _, xsec = pytran.calculate_hitran_xsec(data, M, wns, T=150.,p=10e3)

        print("plotting")
        fig, ax = plt.subplots()
        fig.canvas.set_window_title(molecule)
        ax.semilogy(wns, xsec, 'k-')
        
        ax.set_title(molecule)
        ax.set_ylabel('Cross-section (' + units + ')')
        ax.set_xlabel('wavenumber (cm$^{-1}$)')
        ax.set_xlim([wnmin-0.05*(wnmax-wnmin),wnmax+0.05*(wnmax-wnmin)])


    if False:
        fileAtm='/bira-iasb/projects/planetary/Asimut/CommonData/Atmosphere/Mars/Mars_RefAtmosphere-fd-a461.dat'
        Z, P, T=np.loadtxt(fileAtm, dtype=float, comments='%', usecols=(0,1,2), unpack=True)

        print(len(T), len(wns))

        fig, ax = plt.subplots()

        xsec = np.zeros((len(T),len(wns)))
        for i in range(0,len(T)):
            print("%d of %d" % (i, len(T)))
            xsec[i,:] =  pytran.calculate_hitran_xsec(data, M, wns, T=T[i], p=P[i]*1e3)[1]

            ax.plot(wns, xsec[i,:], label='%.1f'%Z[i])

        ax.set_yscale('log')
        ax.legend()

        with h5py.File('sigma_CH4_134.h5', 'w') as f:
            f['z'] = Z
            f['T'] = T
            f['P'] = P
            f['wvn'] = wns
            f['sigma'] = xsec

    if True:

        fileAtm='/bira-iasb/projects/planetary/Asimut/CommonData/Atmosphere/Mars/Mars_RefAtmosphere-fd-a461.dat'
        Z, P, T=np.loadtxt(fileAtm, dtype=float, comments='%', usecols=(0,1,2), unpack=True)

        Nden = P/T*1e3/(1.380650e-16)
        Ncol = Nden * 10.e5
        print(Nden[-1], Ncol[-1])

        xsec =  pytran.calculate_hitran_xsec(data, M, wns, T=T[-1], p=P[-1]*1e3)[1]

        Trans = np.exp(-xsec*Ncol[-1])

        fig, axs = plt.subplots(2,1)
        axs[0].plot(wns, xsec)
        axs[0].set_yscale('log')
        axs[1].plot(wns, Trans)


    plt.show()
    


    
