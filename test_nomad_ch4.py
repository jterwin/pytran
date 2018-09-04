


## ==================================================================================================
## ==================================================================================================

if (__name__ == "__main__"):

    import os
    import glob
    
    import pytran
    import numpy as np
    from matplotlib import pyplot as plt

    #molecule = 'H2'       ## water
    molecule = 'CO2'       ## carbon dioxide
    #molecule = 'NH3'       ## ammonia
    #molecule = 'SO2'       ## sulfur dioxide
    #molecule = 'CH4'       ## methane
    #molecule = 'H2S'       ## hydrogen sulfide
    #molecule = 'O3'        ## ozone
    #molecule = 'C2H6'      ## ethane

    units = 'cm^2'

    wnmin, wnmax = (3011.44, 3035.44)
    print('(%.1f, %.1f)' % (wnmin,wnmax))

    HITRANDIR = os.path.expanduser('~/lib/HITRAN')
    print(HITRANDIR)

    M = pytran.get_molecule_id(molecule)    
    filenames = glob.glob(os.path.join(HITRANDIR,'HITRAN2012/By-Molecule/Uncompressed-files/%02i_hit*.par' % M))
    filename = filenames[0]
    filename = os.path.normpath(os.path.abspath(filename))
    print(filename)
    
    data = pytran.read_hitran2012_parfile(filename, wnmin, wnmax)
    nlines = len(data['S'])
    print('Found %i lines' % nlines)

    
    ## Next set up for calculating the absorption cross-section, given the HITRAN database's values for the
    ## line centers and line strengths.
    print('Calculating the absorption cross-section spectrum ...')
    (wns, xsec) = pytran.calculate_hitran_xsec(data, M, wnmin, wnmax, T=100.,p=10e1)

     
    print("plotting")
    fig, ax = plt.subplots()
    fig.canvas.set_window_title(molecule)
    ax.semilogy(wns, xsec, 'k-')
    
    ax.set_title(molecule)
    ax.set_ylabel('Cross-section (' + units + ')')
    ax.set_xlabel('wavenumber (cm$^{-1}$)')
    ax.set_xlim([wnmin-0.05*(wnmax-wnmin),wnmax+0.05*(wnmax-wnmin)])




    if False:
        print("sorting")
        fig, ax = plt.subplots()
        ind = xsec.argsort()
        ax.semilogy(np.linspace(0,1,xsec.size),xsec[ind])


        sigmabar = sum(xsec)
        print( '%e' % (sigmabar) )
        sigmabar = sigmabar*(wns[1]-wns[0])  
        print ( '%e' % (sigmabar) )
        print ( '%e' % (sigmabar/(wns[-1]-wns[0])) )
    

    
    plt.show()
    


    
