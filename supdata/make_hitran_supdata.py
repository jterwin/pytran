

from __future__ import print_function


def parse_molparam_file(filename):
    """ read molparam file into dict structure """
    with open(filename) as f:
        lines = f.readlines()
        nmols = 0
        mlist = {}
        for line in lines[1:]:
            vars = line.split()
            if len(vars) == 2:
                nmols += 1
                nisos = 0
                molname = vars[0]
                mlist[nmols] = {}
                mlist[nmols]["name"] = molname
            elif len(vars) == 6:
                nisos += 1
                mlist[nmols][nisos] = {}
                mlist[nmols][nisos]["iso"] = vars[0]
                mlist[nmols][nisos]["abundance"] = vars[1]
                mlist[nmols][nisos]["Q296"] = vars[2]
                mlist[nmols][nisos]["gj"] = vars[3]
                mlist[nmols][nisos]["mass"] = vars[4]
                mlist[nmols][nisos]["gID"] = vars[5]
            else:
                print("skipping line (%d): %s"%(len(vars), line))
    return(mlist)

if __name__ == '__main__':

    import os

    import numpy as np

    testing = False

    mlist = parse_molparam_file('molparam.txt')
    print("read in %d molecules from molparam.txt"%len(mlist))

    NbMol = len(mlist)
    NbIso = np.sum([len(mlist[i+1])-1 for i in range(NbMol)])
    NbMaxIso = np.max([len(mlist[i+1])-1 for i in range(NbMol)])
    print(NbMol, NbIso, NbMaxIso)

    if testing:
        outfile = './hitran_supdate.py'
    else:
        outfile = '../pytran/hitran_supdata.py'

    with open(outfile, 'w') as f:

        f.write("# counts of molecules and isotopologues\n")
        f.write("NbMol = %d\n" % NbMol)
        f.write("NbIso = %d\n" % NbIso)
        f.write("NbMaxIso = %d\n" % NbMaxIso)
        f.write("\n")

        f.write("# molecule parameters\n")
        f.write("molparam = {\n")
        for mol in sorted(mlist):
            f.write("    %d : {\n" % mol)
            f.write("        'name' : '%s',\n" % mlist[mol]['name'])
            for iso in range(20):
                if iso in mlist[mol]:
                    f.write("        %d : {'iso':'%s', 'abundance':%s, 'Q296':%s, 'gj':%s, 'mass':%s},\n" % (iso, mlist[mol][iso]['iso'], mlist[mol][iso]['abundance'], mlist[mol][iso]['Q296'], mlist[mol][iso]['gj'], mlist[mol][iso]['mass']))
            f.write("    },\n")
        f.write("}\n\n")


        f.write("\n")
        f.write("# tabulated qtips values\n")
        f.write("qtab = {\n")
        for mol in sorted(mlist):
            f.write("    %d : {\n" % mol)
            for iso in range(20):
                if iso in mlist[mol]:
                    filename = 'qtips_files/q%s.txt'%mlist[mol][iso]['gID']
                    if os.path.exists(filename):
                        with open(filename, 'r') as f2:
                            lines = f2.readlines()
                        T = np.array([line[:4].strip() for line in lines])
                        q = np.array([line[4:-1].strip() for line in lines])
                        qs = ["%g"%float(val) for val in q[:1000]]
                        #print(mol, iso, T[0], T[-1])
                    else:
                        T = ["%d"%val for val in np.arange(1000) + 1.0]
                        Q296 = float(mlist[mol][iso]['Q296'])
                        qs = ["%g"%Q296 for val in np.ones(1000)] 
                        #print(mol,iso,)
                    
                    #f.write("        %d : {\n" % (iso))
                    #f.write("            'Tmin':%s,\n" % T[0])
                    #f.write("            'Tmax':%s,\n" % '1000.')
                    #f.write("            'q':[%s],\n" % ", ".join(qs))
                    #f.write("        },\n")

                    f.write("        %d : {'Tmin':%s, 'Tmax':%s, 'q':[%s]},\n" % (iso, T[0], '1000', ", ".join(qs)))

            f.write("    },\n")

        f.write("}\n")
        
