
from __future__ import division, print_function, absolute_import

#from numpy import *
import os

__all__ = ['qt_tab', 'molecule_list', 'get_isoname', 'get_ni']

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

    # finish up and return
    for ind in mlist.keys():
        mlist[ind]["mass"] =  sum([a*b for a,b in zip(mlist[ind]["abundances"],mlist[ind]["masses"])])
    return(mlist)


molecule_list = parse_molparam_file(_molparam_file)


def get_isoname(nmol,niso=1):
    """ returns full isotopologue name
    
    Inputs:
        nmol - molecule number
        niso - isotopologue number (default 1)
    
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

def get_ni(name):
    """ return nmol,niso from isotopologue name  """
    
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

#print(get_isoname(6,1))
#print(get_isoname(6,2))
#print(get_isoname(6,3))
#print(get_isoname(6,4))
#print(get_isoname(6,5))
#print(get_isoname(64,1))

#print(get_ni("CO2_838"))
#print(get_ni("CO_38"))
#print(get_ni("CO2_222"))
#print(get_ni("C6H12O6_838"))



full_iso_list = []
for nmol in range(1,len(molecule_list)+1):
    for niso in range(len(molecule_list[str(nmol)]["isops"])):
        full_iso_list.append(get_isoname(nmol,niso+1))

#print(full_iso_list)
#print(set(full_iso_list)-set(_qt_tab_2003.keys()))
#print(set(_qt_tab_2003.keys())-set(full_iso_list))

for i in range(len(full_iso_list)):
    name = full_iso_list[i]
    #print(name, _qt_tab_2003[name])
    #print(name, qt_tab[name])

