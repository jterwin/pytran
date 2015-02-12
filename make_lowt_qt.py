from __future__ import division, print_function, absolute_import

from numpy import *
import scipy.optimize as optimization
import os

from hitran_qt import *


full_iso_list = []
for nmol in range(1,len(molecule_list)+1):
    for niso in range(len(molecule_list[str(nmol)]["isops"])):
        full_iso_list.append(get_isoname(nmol,niso+1))


def rotq_spherical_top(tmp, B0 = 5.2410356, D0 = 1.10864e-4, Ig = 0.5, csn = 12., n = 4):
    """  rotational partition function for spherical top molecule
    adapted from:
      McDowell, Robin S. Rotational Partition Function for Spherical-Top Molecules. Journal
      Quant. Spectrosc. Radiat. Transfer, Vol 38, 1987.
    """

    #from numpy import exp, pi
    import scipy.constants as constants
    h = constants.physical_constants['Planck constant'][0]*1e7
    kb = constants.physical_constants['Boltzmann constant'][0]*1e7
    c = constants.physical_constants['speed of light in vacuum'][0]*1e2
    amu = constants.physical_constants['atomic mass constant'][0]*1e3
                       
    sigma_star = (2*Ig+1)**n/csn
    beta = h*c*B0/kb/tmp
    delta = ((16.*pi/3.**1.5)*exp(-pi**2/9./beta)*(1.+2.*exp(-pi**2/3./beta)) +
             3.*pi*exp(-pi**2/4./beta))/(2.*Ig+1.)**2
    Qr = sigma_star*pi**0.5*exp(beta/4)*beta**-1.5*(1.+delta)*(1+15.*D0/(4.*beta*B0))

    return(Qr)

def rotq_symmetric_top(tmp, A = 5.250774, B = 3.880159, Ig = 0.5, csn = 3., n = 3, p = 1./2.,
                       rho0 = -1.191e-5, rho1 = 5.97683e-5, rho2 = 0.46345e-9, rho3 = -1.415e-12):
    """  rotational partition function for symmetric top molecule
    adapted from:
      McDowell, Robin S. Rotational Partition Function for Symmetric-Top Molecules. The Journal
      of Chemical Physics, Vol 93, 1990.
      
    A - centrifugal distortion constant (cm^-1)
    B - rotational constant (cm^-1)
    Ig - nuclear spin number
    csn - classical symmetry number
    n - number of X atom in WX_3Y molecule
    p = 2/(2I+1)^2, where I is spin of X atom in WX_3Y molecule
    rho0,rho1,rho2,rho3 - polynomial correction for centrifugal distortion
    """

    import scipy.constants as constants
    h = constants.physical_constants['Planck constant'][0]*1e7
    kb = constants.physical_constants['Boltzmann constant'][0]*1e7
    c = constants.physical_constants['speed of light in vacuum'][0]*1e2
    amu = constants.physical_constants['atomic mass constant'][0]*1e3
 
    m = B/A
    beta = h*c*B/kb/tmp
    sigma_star = (2*Ig+1)**n/csn

    delta = p*exp(-pi**2*m/n**2/beta)*(1. + exp(pi**2*(2.-n)*m/n/beta) +
                pi**2*m/90./n**2*(15.+4*beta*(1-m)+15.*(1-n)**2*exp(pi**2*(2-n)*m/n/beta)) +
                7.*pi**4*m**4/360./n**4    )
    Qr = sigma_star*((pi*m)**0.5*exp(beta*(4-m)/12.)*beta**-1.5*(1.+beta**2*(1-m)**2/90.)*
                     (1.+delta)*(1. + rho0 + rho1*beta**-1 + rho2*beta**-2 + rho3*beta**-3))
    
    return(Qr)



def rotq_linear(tmp, B = 1.478221825, D = 2.909851e-6, H = 0., I = 1., csn = 1., kappa = 0.):
    """ rotational partition function for linear molecules
    adapted from:
        McDowell, Robin S. Rotational partition functions for linear molecules

    B - rotational constant (cm^-1)
    D - centrifugal distortion quartic constant (cm^-1)
    H - centrifugal distortion sextic constant (cm^-1)
    I - nuclear spin multiplicity (may be compensated by hitran molecule multiplicity)
    csn - classical symmetic number
    kappa - nuclear spin correction
    """

    import scipy.constants as constants
    h = constants.physical_constants['Planck constant'][0]*1e7
    kb = constants.physical_constants['Boltzmann constant'][0]*1e7
    c = constants.physical_constants['speed of light in vacuum'][0]*1e2
    amu = constants.physical_constants['atomic mass constant'][0]*1e3
    
    dd = D/B
    hh = H/B
    beta = h*c*B/kb/tmp

    fc = 1. + 2.*dd*(3-beta)/(3.*beta) + 6.*(2.*dd**2-hh)/(beta**2) + 120.*dd*(dd**2-hh)/(beta**3)
    Qr = csn**-1*I**2*exp(beta/3.)*beta**-1*(1. + beta**2/90. + 8.*beta**3/2835. +
                                             kappa*I**-1*pi**1.5*exp(-beta/12.-pi**2/(4.*beta))*beta**-0.5)*fc
    
    return(Qr)


lowtlist = []
qt_lowT = {}
tmps = [10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,85.0]

qt_lowT["Temp"] = tmps
for i in range(len(full_iso_list)):
    name = full_iso_list[i]

    if name in ("CH4_211"):
        #print(name)
        lowtlist.append(name)
        qt_lowT[name] = []
        for tmp in tmps:
            qt_lowT[name].append(rotq_spherical_top(tmp, B0 = 5.2410356, D0 = 1.10864e-4,
                                                    Ig = 0.5, csn = 12., n = 4))
        #print(qt_lowT[name])

    if name in ("CH4_311"):
        #print(name)
        lowtlist.append(name)
        qt_lowT[name] = []
        for tmp in tmps:
            qt_lowT[name].append(rotq_spherical_top(tmp, B0 = 5.241253, D0 = 1.1072e-4,
                                                    Ig = 0.5, csn = 12., n = 4))
        
        nmol,niso = get_ni(name)
        gj = molecule_list[str(nmol)]["gj"][niso-1]
        qt_lowT[name] = gj*array(qt_lowT[name])
        
    if name in ("CH4_212","CH4_312"):
        #print(name)
        lowtlist.append(name)
        qt_lowT[name] = []
        for tmp in tmps:
            qt_lowT[name].append(
                rotq_symmetric_top(tmp, A = 5.250774, B = 3.880159, Ig = 0.5, csn = 3., n = 3, p = 1./2.,
                       rho0 = -1.191e-5, rho1 = 5.97683e-5, rho2 = 0.46345e-9, rho3 = -1.415e-12) )
        
        nmol,niso = get_ni(name)
        gj = molecule_list[str(nmol)]["gj"][niso-1]
        qt_lowT[name] = gj*array(qt_lowT[name])
        

    if name in ("HCN_124"):
        print(name)
        lowtlist.append(name)
        
        qt_lowT[name] = []
        for tmp in tmps:
            qt_lowT[name].append(
                rotq_linear(tmp, B = 1.478221825, D = 2.909851e-6, H = 2.72e-12, I = 1., csn = 1., kappa = 0.) )
        
        nmol,niso = get_ni(name)
        gj = molecule_list[str(nmol)]["gj"][niso-1]
        qt_lowT[name] = gj*array(qt_lowT[name])

    if name in ("HCN_134"):
        print(name)
        lowtlist.append(name)
        
        qt_lowT[name] = []
        for tmp in tmps:
            qt_lowT[name].append(
                rotq_linear(tmp, B = 1.439999, D = 2.769e-6, H = 0., I = 1., csn = 1., kappa = 0.) )
        
        nmol,niso = get_ni(name)
        gj = molecule_list[str(nmol)]["gj"][niso-1]
        qt_lowT[name] = gj*array(qt_lowT[name])

    if name in ("HCN_125"):
        print(name)
        lowtlist.append(name)
        
        qt_lowT[name] = []
        for tmp in tmps:
            qt_lowT[name].append(
                rotq_linear(tmp, B = 1.435249, D = 2.745e-6, H = 0., I = 1., csn = 1., kappa = 0.) )
        
        nmol,niso = get_ni(name)
        gj = molecule_list[str(nmol)]["gj"][niso-1]
        qt_lowT[name] = gj*array(qt_lowT[name])

    if name in ("CO_26"):
        print(name)
        lowtlist.append(name)
        
        qt_lowT[name] = []
        for tmp in tmps:
            qt_lowT[name].append(
                rotq_linear(tmp, B = 1.92252896, D = 6.121095e-6, H = 5.73e-12, I = 1., csn = 1., kappa = 0.) )
        
        nmol,niso = get_ni(name)
        gj = molecule_list[str(nmol)]["gj"][niso-1]
        qt_lowT[name] = gj*array(qt_lowT[name])

    if name in ("CO_36"):
        print(name)
        lowtlist.append(name)
        
        qt_lowT[name] = []
        for tmp in tmps:
            qt_lowT[name].append(
                rotq_linear(tmp, B = 1.83797198, D = 5.59347e-6, H = 5.06e-12, I = 1., csn = 1., kappa = 0.) )
        
        nmol,niso = get_ni(name)
        gj = molecule_list[str(nmol)]["gj"][niso-1]
        qt_lowT[name] = gj*array(qt_lowT[name])

    if name in ("CO_28"):
        print(name)
        lowtlist.append(name)
        
        qt_lowT[name] = []
        for tmp in tmps:
            qt_lowT[name].append(
                rotq_linear(tmp, B = 1.83098081, D = 5.55071e-6, H = 4.44e-12, I = 1., csn = 1., kappa = 0.) )
        
        nmol,niso = get_ni(name)
        gj = molecule_list[str(nmol)]["gj"][niso-1]
        qt_lowT[name] = gj*array(qt_lowT[name])

    if name in ("CO_27"):
        print(name)
        lowtlist.append(name)
        
        qt_lowT[name] = []
        for tmp in tmps:
            qt_lowT[name].append(
                rotq_linear(tmp, B = 1.87396018, D = 5.81502e-6, H = 0., I = 1., csn = 1., kappa = 0.) )
        
        nmol,niso = get_ni(name)
        gj = molecule_list[str(nmol)]["gj"][niso-1]
        qt_lowT[name] = gj*array(qt_lowT[name])

    if name in ("CO_38"):
        print(name)
        lowtlist.append(name)
        
        qt_lowT[name] = []
        for tmp in tmps:
            qt_lowT[name].append(
                rotq_linear(tmp, B = 1.7464085, D = 5.04866e-6, H = 0., I = 1., csn = 1., kappa = 0.) )
        
        nmol,niso = get_ni(name)
        gj = molecule_list[str(nmol)]["gj"][niso-1]
        qt_lowT[name] = gj*array(qt_lowT[name])

    if name in ("CO_37"):
        print(name)
        lowtlist.append(name)
        
        qt_lowT[name] = []
        for tmp in tmps:
            qt_lowT[name].append(
                rotq_linear(tmp, B = 1.7922, D = 5.3211e-6, H = 0., I = 1., csn = 1., kappa = 0.) )
        
        nmol,niso = get_ni(name)
        gj = molecule_list[str(nmol)]["gj"][niso-1]
        qt_lowT[name] = gj*array(qt_lowT[name])
        

with open('data/parsum_lowT.dat', 'w') as f:
    print("")
    line = ["Temp"]
    for name in lowtlist:
        line.append(name)
    print("".join([t.rjust(20) for t in line]))
    f.write("".join([t.rjust(20) for t in line]))
    
    for i in range(len(qt_lowT["Temp"])):
        line = []
        line.append("%20.1f"%qt_lowT["Temp"][i])
        for name in lowtlist:
            line.append("%20.5e"%qt_lowT[name][i])
        print("".join(line))
        f.write("\n"+"".join(line))
    f.write("\n")
        
        

print("")
for i in range(5):
    line = []
    line.append("%20.1f"%qt_tab["Temp"][i])
    for name in lowtlist:
        line.append("%20.5e"%qt_tab[name][i])
    print("".join(line))




def print_array(array, fmt=None):

    if fmt is None:
        fmt = "%10f"
    print("".join([fmt%x for x in array]))


#
# check fitting routine for ch4_211
# in the end, i recommend method 2 for spherical symmetric molecules
#
if (False):
    name = "CH4_211"
    print(name)
    xdata, ydata = qt_tab["Temp"][1:5], qt_tab[name][1:5]
    print(xdata,ydata)

    def func1(params, x, y):
        """ sigma_star = params[0], theta = params[1] """
        beta = params[1]/x
        return (y-params[0]*pi**0.5*exp(beta/4)*beta**-1.5)
    def func2(params, x, y):
        """ sigma_star = params[0], theta = params[1] """
        beta = params[1]/x
        delta = (4.*pi/3.**1.5)*exp(-pi**2/9./beta)
        return (y-params[0]*pi**0.5*exp(beta/4)*beta**-1.5*(1.+delta))
    def func3(params, x, y):
        """ sigma_star = params[0], theta = params[1], m = D/B = params[2]"""
        beta = params[1]/x
        delta = (4.*pi/3.**1.5)*exp(-pi**2/9./beta)
        return (y-params[0]*pi**0.5*exp(beta/4)*beta**-1.5*(1.+delta)*(1.+15./4./beta*params[2]))

    sol1,J = optimization.leastsq(func1, (10., 1.), args=(xdata,ydata))
    print(sol1,J)
    sol2,J = optimization.leastsq(func2, (10., 1.), args=(xdata,ydata))
    print(sol2,J)
    sol3,J = optimization.leastsq(func3, (10., 1., 0.), args=(xdata,ydata))
    print(sol3,J)

    print("")
    print_array(xdata,fmt="%14.3f")
    print_array(ydata,fmt="%14.3f")
    beta = sol1[1]/xdata
    delta = (4.*pi/3.**1.5)*exp(-pi**2/9./beta)
    print_array(sol1[0]*pi**0.5*exp(beta/4)*beta**-1.5*(1.+delta),fmt="%14.3f")

    print("")
    print_array(qt_lowT["Temp"],fmt="%14.3f")
    print_array(qt_lowT[name],fmt="%14.3f")
    beta = sol1[1]/qt_lowT["Temp"]
    print_array(sol1[0]*pi**0.5*exp(beta/4)*beta**-1.5,fmt="%14.3f")
    beta = sol2[1]/qt_lowT["Temp"]
    delta = (4.*pi/3.**1.5)*exp(-pi**2/9./beta)
    print_array(sol2[0]*pi**0.5*exp(beta/4)*beta**-1.5*(1.+delta),fmt="%14.3f")
    beta = sol3[1]/qt_lowT["Temp"]
    delta = (4.*pi/3.**1.5)*exp(-pi**2/9./beta)
    print_array(sol3[0]*pi**0.5*exp(beta/4)*beta**-1.5*(1.+delta)*(1.+15./4./beta*sol3[2]),fmt="%14.3f")


if (False):
    name = "HCN_124"
    print(name)
    xdata, ydata = qt_tab["Temp"][1:5], qt_tab[name][1:5]
    print(xdata,ydata)

    def func1(params, x, y):
        """ sigma_star = params[0], theta = params[1] """
        beta = params[1]/x
        return (y-params[0]*exp(beta/3.)*beta**-1.0)
    #def func2(params, x, y):
    #    """ sigma_star = params[0], theta = params[1] """
    #    beta = params[1]/x
    #    delta = (4.*pi/3.**1.5)*exp(-pi**2/9./beta)
    #    return (y-params[0]*pi**0.5*exp(beta/4)*beta**-1.5*(1.+delta))
    #def func3(params, x, y):
    #    """ sigma_star = params[0], theta = params[1], m = D/B = params[2]"""
    #    beta = params[1]/x
    #    delta = (4.*pi/3.**1.5)*exp(-pi**2/9./beta)
    #    return (y-params[0]*pi**0.5*exp(beta/4)*beta**-1.5*(1.+delta)*(1.+15./4./beta*params[2]))

    sol1,J = optimization.leastsq(func1, (10., 1.), args=(xdata,ydata))
    print(sol1,J)
    #sol2,J = optimization.leastsq(func2, (10., 1.), args=(xdata,ydata))
    #print(sol2,J)
    #sol3,J = optimization.leastsq(func3, (10., 1., 0.), args=(xdata,ydata))
    #print(sol3,J)

    print("")
    print_array(xdata,fmt="%14.3f")
    print_array(ydata,fmt="%14.3f")
    beta = sol1[1]/xdata
    print_array(sol1[0]*exp(beta/3.)*beta**-1.0,fmt="%14.3f")

    print("")
    print_array(qt_lowT["Temp"],fmt="%14.3f")
    print_array(qt_lowT[name],fmt="%14.3f")
    beta = sol1[1]/qt_lowT["Temp"]
    print_array(sol1[0]*exp(beta/3.)*beta**-1.0,fmt="%14.3f")
    #beta = sol2[1]/qt_lowT["Temp"]
    #delta = (4.*pi/3.**1.5)*exp(-pi**2/9./beta)
    #print_array(sol2[0]*pi**0.5*exp(beta/4)*beta**-1.5*(1.+delta),fmt="%14.3f")
    #beta = sol3[1]/qt_lowT["Temp"]
    #delta = (4.*pi/3.**1.5)*exp(-pi**2/9./beta)
    #print_array(sol3[0]*pi**0.5*exp(beta/4)*beta**-1.5*(1.+delta)*(1.+15./4./beta*sol3[2]),fmt="%14.3f")

