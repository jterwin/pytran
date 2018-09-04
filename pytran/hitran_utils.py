
from __future__ import division, print_function, absolute_import

from pytran.hitran_supdata import molparam, qtab

__all__ = ['get_molecule_id', 'get_iso_id', 'get_molecule_mass', 'get_iso_name',
           'get_iso_mass', 'get_molecule_nisops', 'qtips']


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
    
    for mol in molparam:
        if (molparam[mol]["name"] == name):
            nmol = mol
            break

    if (nmol == 0):
        raise LookupError("get_ni: molecule %s not found"%name)
    
    return nmol
    
    
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

    for mol in molparam:
        if (molparam[mol]["name"] == smol):
            nmol = int(mol)
            break

    if (nmol == 0):
        raise LookupError("get_ni: molecule %s not found"%smol)
    
    ind = str(nmol)
    for i in range(1,len(molparam[nmol])):
        if (molparam[nmol][i]['iso'] == siso):
            niso = i
            break

    if (niso == 0):
        raise LookupError("get_ni: isotopologue %s of %s not found"%(siso,smol))
    
    return (nmol,niso)


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
    
    if (nmol in molparam):
        if (niso in molparam[nmol]):
            isoname = "%s_%s" % (molparam[nmol]["name"], molparam[nmol][niso]['iso'])
        else:
            raise LookupError("ni_to_isoname: niso=%d for %s not found"%
                              (niso,molparam[nmol]["name"]))
    else:
        #isoname = None
        raise LookupError("ni_to_isoname: nmol=%d is not found"%(nmol))

    return isoname


def get_molecule_nisops(nmol):
    """
    Returns number of isotopologues for molecule

    Parameters
    ----------
    nmol : integer
        HITRAN molecule number

    Returns
    -------
    nisop : int
        number of isotopologues of molecule

    """
    
    if nmol not in molparam:
        raise LookupError("nmol=%d not in molparam"%nmol)

    return (len(molparam[nmol])-1)

        
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

    if nmol not in molparam:
        raise LookupError("nmol=%d not in molparam"%nmol)

    nisops = get_molecule_nisops(nmol)
    mass = sum([molparam[nmol][niso]["abundance"]*molparam[nmol][niso]["mass"] for niso in range(1,nisops+1)])

    return mass


def get_iso_mass(nmol,niso=1):
    """
    Returns  mass for isotopologue

    Parameters
    ----------
    nmol : integer
        HITRAN molecule number
    niso : integer
        isotopologue number

    Returns
    -------
    mass : float
        mass of isotopologue

    """

    if nmol not in molparam:
        raise LookupError("nmol=%d not in molparam"%nmol)
    if niso not in molparam[nmol]:
        raise LookupError("niso=%d not in molparam[%d]"%(niso,nmol))

    mass = molparam[nmol][niso]["mass"]
    return mass
    

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

    if nmol not in qtab:
        raise LookupError("nmol=%d not in qtab"%nmol)
    if niso not in qtab[nmol]:
        raise LookupError("niso=%d not in qtab[%d]"%(niso,nmol))
    if tmp < qtab[nmol][niso]['Tmin'] or tmp > qtab[nmol][niso]['Tmax']:
        raise ValueError("tmp=%f outside of range (%.1f, %.1f)" % (tmp, qtab[nmol][niso]['Tmin'], qtab[nmol][niso]['Tmax']))

    dt = tmp - qtab[nmol][niso]['Tmin']
    it = int(dt)
    ft = dt - it
    q = qtab[nmol][niso]['q'][it]*(1.-ft) + qtab[nmol][niso]['q'][it+1]*ft

    return q
