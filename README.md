# pytran

A python package to read and use the HITRAN database.

### Intro
This is my attempt to create a simple interface into the HITRAN database for the purpose of using IR absorption cross sections. 

One should obtain their own copy of hitran parameter files, for instance from the hitran website (https://hitran.org/) or using hapy (https://hitran.org/hapi/).

### Usage

#### reading linelist and computing cross section



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

