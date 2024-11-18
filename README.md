# pytran

A python package to read and use the HITRAN database.

### Intro
This is my attempt to create a simple interface into the HITRAN database for the purpose of using IR absorption cross sections. 

One should obtain their own copy of hitran parameter files, for instance from the hitran website (https://hitran.org/) or using hapy (https://hitran.org/hapi/).

### Usage

#### reading linelist and computing cross section

In the simplest example, you can read in a par file to a linelist dictionary (using `pytran.read_hitran2012_parfile`) and then compute the extiction cross section (using `pytran.calculate_hitran_xsec`):
```
    parfilename = './supdata/05_hit20_0_1000.par'
    linelist = pytran.read_hitran2012_parfile(parfilename)

    wn = np.arange(0., 400., 1.e-3)
    sigma = pytran.calculate_hitran_xsec(linelist, wn, T, P)
```
The script `test_hitran_lbl.py` provides a complete example.

Both function can take additional arguments refine the computation. 

#### utilities

Several functions are provided which interact with the molecular parameters and total internal partition sums (`get_molecule_id`, `get_iso_id`, `get_molecule_mass`, `get_iso_name`, `get_iso_mass`, `get_molecule_nisops`, `qtips`). They are mainly used by `pytran.calculate_hitran_xsec` to correct the line strength and line parameters. 

The script `test_supdata.py` provides some examples on their usages.

### Design Goals

##### hitran.py
- [x] Read HITRAN .par files into linelist dictionary
- [x] Compute cross section from linelist dictionary, correcting for temperature and pressure conditions
- [x] Use partition function in computing cross section (updated to 2020)
- [x] Use isotopologue values to compute line-strength and spectra
- [ ] Create new function 'read_linelist', have it switch between different filetypes (e.g. par file, hdf5)
- [ ] Add utility to 'read_linelist' to split by isotopologue and/or split by upper/lower global quanta (return dict, or nested dict, of linelist)


##### hitran_supplemental.py
- [x] Read HITRAN molparam.txt file into molecule_param dictionary
- [x] Create utilities to return names, indicies, values for molparam
dict
- [x] Read HITRAN parsum.dat into parsum dictionary
- [x] Create utility to evaluate parsum for particular isotopologue

##### todo
- [x] Develop into a Python Package
- [ ] Add line-mixing for CO2 and CH4
- [ ] Add hitran CIA
- [ ] Create utilities to interpret linelists (i.e. list of transitions with bounds and strengths, compute band Einstein coefficient, etc)
- [ ] Create utility to create k-coefficient table
- [ ] Make a GUI like pyHAWKS (probably not)

### Aknowledgements

I adapted the initial code from https://github.com/nzhagen/hitran, thanks go to the author.

A sample line list for CO rotational lines is provided for the example. The references are also included. 

