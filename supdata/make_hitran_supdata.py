

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
                #mlist[ind]["isops"] = []
                #mlist[ind]["abundances"] = []
                #mlist[ind]["Q296"] = []
                #mlist[ind]["gj"] = []
                #mlist[ind]["masses"] = []
            elif len(vars) == 6:
                nisos += 1
                mlist[nmols][nisos] = {}
                

                mlist[nmols][nisos]["iso"] = vars[0]
                mlist[nmols][nisos]["abundance"] = vars[1]
                mlist[nmols][nisos]["Q296"] = vars[2]
                mlist[nmols][nisos]["gj"] = vars[3]
                mlist[nmols][nisos]["mass"] = vars[4]
                mlist[nmols][nisos]["gID"] = vars[5]
    return(mlist)

if __name__ == '__main__':

    import os

    import numpy as np
    from matplotlib import pyplot as plt    


    mlist = parse_molparam_file('molparam.txt')

    with open('../pytran/hitran_supdata.py', 'w') as f:

        f.write("molparam = {}\n")
        for mol in sorted(mlist):
            f.write("molparam[%d] = {\n" % mol)
            f.write("    'name' : '%s',\n" % mlist[mol]['name'])
            for iso in range(20):
                if iso in mlist[mol]:
                    f.write("    %d : {'iso':'%s', 'abundance':%s, 'Q296':%s, 'gj':%s, 'mass':%s},\n" % (iso, mlist[mol][iso]['iso'], mlist[mol][iso]['abundance'], mlist[mol][iso]['Q296'], mlist[mol][iso]['gj'], mlist[mol][iso]['mass']))
            f.write("}\n")


        f.write("qtab = {}\n")
        for mol in sorted(mlist):
            f.write("qtab[%d] = {}\n" % mol)
            for iso in range(20):
                if iso in mlist[mol]:
                    filename = 'qtips_files/q%s.txt'%mlist[mol][iso]['gID']
                    if os.path.exists(filename):
                        with open(filename, 'r') as f2:
                            lines = f2.readlines()
                        T = np.array([line[:4].strip() for line in lines])
                        q = np.array([line[4:-1].strip() for line in lines])
                    else:
                        T = ["%d"%val for val in np.arange(1000) + 1.0]
                        q = ["%d"%val for val in np.ones(1000)]
                        print(mol,iso,len(T),len(q))
                    
                    f.write("qtab[%d][%d] = {\n" % (mol, iso))
                    f.write("    'Tmin':%s,\n" % T[0])
                    f.write("    'Tmax':%s,\n" % T[-1])
                    f.write("    'q':[%s],\n" % ", ".join(q))
                    f.write("}\n")

    #fig, ax = plt.subplots()    

    for i in range(1):
        
        filename = 'qtips_files/q%d.txt'%i
        #print(filename)

        if os.path.exists(filename):
            #T, q = np.loadtxt(filename, unpack=True)

            #T = np.zeros(10000)
            #q = np.zeros(10000)
            with open(filename, 'r') as f:
                lines = f.readlines()
                T = np.array([float(line[:4]) for line in lines])
                q = np.array([float(line[4:]) for line in lines])
                
        else:
            continue

        print(i, len(T), T[0], T[-1], q[0], q[-1])

        #ax.plot(T/(np.arange(len(T))+T[0]))

    #plt.show()

        
