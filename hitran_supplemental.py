
from __future__ import division, print_function, absolute_import

#from numpy import *
import os

__all__ = ['get_molecule_id', 'get_iso_id', 'get_molecule_mass', 'get_iso_name', 'qtips']

_parsum_file_2000 = os.path.expanduser("data/parsum.dat")
_parsum_file_2003 = os.path.expanduser("data/parsum2003.dat")
_parsum_file_2011 = os.path.expanduser("data/parsum2011.dat")

_molparam_file = os.path.expanduser("data/molparam.txt")


def parse_qt_file(filename):
    """ read parsum file into dict structure """
    with open(filename) as f:
        lines = f.readlines()
        
        # create dict fields
        tab={}
        colnames = lines[0].split()
        for name in colnames:
            tab[name]=[]
            
        #colnames = lines[0].split()
        #tab = dict.fromkeys(colnames,[])

        # fill in values
        for i in range(1,len(lines)):
            vals = map(float, lines[i].split())
            for j,name in enumerate(colnames):
                tab[name].append(vals[j])
    
    # rename temp field
    tab["Temp"] = tab["Temp(K)"]
    del tab["Temp(K)"]
    
    return(tab)




def parse_molparam_file(filename):
    """ read molparam file into dict structure """

    with open(filename) as f:
        lines = f.readlines()

        nmols = 0
        mlist = {}
        for line in lines[1:]:
            vars = line.split()
        
            if len(vars)==2:
                nmols = nmols+1
                ind = str(nmols)
                mlist[ind] = {}
                mlist[ind]["name"] = vars[0]
                mlist[ind]["isops"] = []
                mlist[ind]["abundances"] = []
                mlist[ind]["Q296"] = []
                mlist[ind]["gj"] = []
                mlist[ind]["masses"] = []
                
            elif len(vars)==5:
                mlist[ind]["isops"].append(vars[0])
                mlist[ind]["abundances"].append(float(vars[1]))
                mlist[ind]["Q296"].append(float(vars[2]))
                mlist[ind]["gj"].append(int(vars[3]))
                mlist[ind]["masses"].append(float(vars[4]))
    
    return(mlist)





#_qt_tab_2000 = parse_qt_file(_parsum_file_2000)
_qt_tab_2003 = parse_qt_file(_parsum_file_2003)
#_qt_tab_2011 = parse_qt_file(_parsum_file_2011)


#print(set(_qt_tab_2003.keys())-set(_qt_tab_2000.keys()))
#print(set(_qt_tab_2011.keys())-set(_qt_tab_2000.keys()))
#print(set(_qt_tab_2011.keys())-set(_qt_tab_2003.keys()))
#print(set(_qt_tab_2003.keys())-set(_qt_tab_2011.keys()))

# set qt_tab
qt_tab = {}
qt_tab.update(_qt_tab_2003)
#qt_tab.update(_qt_tab_2011)  # unecessary

#print(set(qt_tab.keys())-set(_qt_tab_2003.keys()))
#print(set(qt_tab.keys())-set(_qt_tab_2011.keys()))
#print(set(_qt_tab_2003.keys())-set(qt_tab.keys()))
#print(set(_qt_tab_2011.keys())-set(qt_tab.keys()))


molecule_list = parse_molparam_file(_molparam_file)




def get_molecule_id(name):
    """
    For a given input molecular formula, return the corresponding HITRAN
    molecule identifier number.

    Parameters
    ----------
    name : str
        The string describing the molecule.

    Returns
    -------
    M : int
        The HITRAN molecular identified number.
        
    """

    nmol = 0
    
    for ind in molecule_list.keys():
        if (molecule_list[ind]["name"] == name):
            nmol = int(ind)
            break

    if (nmol == 0):
        raise LookupError("get_ni: molecule %s not found"%name)
    
    return(nmol)
    
    
def get_iso_id(name):
    """ 
    For a given input molecular formula, return the corresponding HITRAN molecule
    and isotopologue identifier numbers.

    Parameters
    ----------
    name : str
        The string describing the molecule.

    Returns
    -------
    nmol : int
        The HITRAN molecular identifier number.
    niso : int
        The HITRAN isotopologue identifier number.
        
    """
    
    nmol, niso = 0, 0
    smol, siso = name.split('_')

    for ind in molecule_list.keys():
        if (molecule_list[ind]["name"] == smol):
            nmol = int(ind)
            break

    if (nmol == 0):
        raise LookupError("get_ni: molecule %s not found"%smol)
    
    ind = str(nmol)
    for i in range(len(molecule_list[ind]["isops"])):
        if (molecule_list[ind]["isops"][i] == siso):
            niso = i+1
            break

    if (niso == 0):
        raise LookupError("get_ni: isotopologue %s of %s not found"%(siso,smol))
    
    return(nmol,niso)

    
def get_molecule_mass(nmol):
    """
    Returns average mass for molecule

    Parameters
    ----------
    nmol : integer
        HITRAN molecule number

    Returns
    -------
    mass : float
        average mass of molecule

    """
    mass = sum([a*b for a,b in zip(molecule_list[str(nmol)]["abundances"],molecule_list[str(nmol)]["masses"])])
    return(mass)

    
def get_iso_name(nmol,niso=1):
    """
    Returns full isotopologue name
    
    Parameters
    ----------
    nmol : integer
        molecule number
    niso : integer
        isotopologue number (default 1)

    Returns
    -------
    isoname : string
        full isotopologue name

    Notes
    -----
        Rasises LookupError Exception if nmol doesn't exist, or niso is out of range for niso, return None
    
    """
    
    if (str(nmol) in molecule_list):
        try:
            isoname = "%s_%s" % (molecule_list[str(nmol)]["name"],
                                 molecule_list[str(nmol)]["isops"][niso-1])
        except IndexError:
            #isoname = None
            raise LookupError("ni_to_isoname: niso=%d for %s not found"%
                              (niso,molecule_list[str(nmol)]["name"]))
    else:
        #isoname = None
        raise LookupError("ni_to_isoname: nmol=%d is not found"%(nmol))
    return(isoname)



def polint4(xx,yy,x):
    """
    Use four point lagrange interpolation to find a value at x
    
    Parameters
    ----------
    xx : float array
        array of x values
    yy : float array
        array of y values
    x : float 
        target x value

    Returns
    -------
    y : float
        interpolate y value
    
    """
    
    from bisect import bisect_left
    from scipy.interpolate import lagrange

    if (x <= xx[0]):
        y = yy[0]
    elif (x >= xx[-1]):
        y = yy[-1]
    else:
        i = bisect_left(xx,x)
        nt1 = min(max(i-2, 0), len(xx)-4)
        nt2 = nt1+4
        c = lagrange(xx[nt1:nt2],yy[nt1:nt2])
        y = c(x)

    return y


def qtips(tmp, nmol, niso=1):
    """
    Computes TIPS value HITRAN molecule
    
    Parameters
    ----------
    tmp : float
        temperature
    nmol : integer
        molecule number
    niso : integer
        isotopologue number (default 1)

    Returns
    -------
    q : float
        interpolated TIPS value

    Notes
    -----
    Uses parsum values for 

    """

    isoname = get_iso_name(nmol,niso)

    if (tmp < 20.0):
        raise ValueError("qtips: tmp < 20 not supported")
    elif (tmp < 70.0):
        raise ValueError("qtips: tmp < 70 will be supported soon")
    elif (tmp < 3000.0):
        print(tmp, nmol, niso, isoname)
        q = polint4(qt_tab["Temp"], qt_tab[isoname], tmp)
    else:
        raise ValueError("qtips: tmp > 3000 not supported")

    return q
    





if (__name__ == "__main__"):

    #print(molecule_list['1'])

    '''
    print(get_molecule_id('CO2'))
    print(get_molecule_id('CH4'))
    print(get_molecule_id('CO'))
    print(get_molecule_id('HCN'))
    print(get_molecule_id('H2O'))
    '''

    '''
    print(get_iso_name(6,1))
    print(get_iso_name(6,2))
    print(get_iso_name(6,3))
    print(get_iso_name(6,4))
    #print(get_iso_name(6,5))
    #print(get_iso_name(64,1))
    '''

    '''
    print(get_iso_id("CO2_838"))
    print(get_iso_id("CO_38"))
    #print(get_iso_id("CO2_222"))
    #print(get_iso_id("C6H12O6_838"))
    '''


    print(qtips(100,6,1))


    full_iso_list = []
    for nmol in range(1,len(molecule_list)+1):
        for niso in range(len(molecule_list[str(nmol)]["isops"])):
            full_iso_list.append(get_iso_name(nmol,niso+1))

    #print(full_iso_list)
    #print(set(full_iso_list)-set(_qt_tab_2003.keys()))
    #print(set(_qt_tab_2003.keys())-set(full_iso_list))

    for i in range(len(full_iso_list)):
        name = full_iso_list[i]
        #print(name, _qt_tab_2003[name])
        #print(name, qt_tab[name])

    
