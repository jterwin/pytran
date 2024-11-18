

if (__name__ == "__main__"):

    import pytran
    import numpy as np
    from matplotlib import pyplot as plt

    parfilename = './supdata/05_hit20_0_1000.par'
    linelist = pytran.read_hitran2012_parfile(parfilename)


    wn = np.arange(0., 400., 1.e-3)
    sigma = pytran.calculate_hitran_xsec(linelist, wn)

    fig, ax = plt.subplots()
    ax.plot(wn, sigma)
    ax.set_yscale('log')

    ax.set_xlabel(r'$\nu$ (cm$^{-1}$)')
    ax.set_ylabel(r'(c$^2$)')


    wn2 = np.arange(20.,50.,1e-3)
    Nbwn = len(wn2)
    Nbiso = pytran.get_molecule_nisops(5)
    sigma2 = np.zeros((Nbiso,Nbwn))

    for iso in range(Nbiso):
        linelist = pytran.read_hitran2012_parfile(parfilename, isolist=(iso+1,), Smin=1.e-28)
        sigma2[iso,:] = pytran.calculate_hitran_xsec(linelist, wn2, T=150., P=7., Pref=1013.)

    fig, ax = plt.subplots()
    ax.plot(wn2, sigma2.sum(axis=0), 'k', alpha=0.5)
    for iso in range(Nbiso):
        ax.plot(wn2, sigma2[iso,:], label=pytran.get_iso_name(5,iso+1))
    ax.set_yscale('log')
    ax.legend()
    ax.set_xlabel(r'$\nu$ (cm$^{-1}$)')
    ax.set_ylabel(r'(c$^2$)')

    plt.show()