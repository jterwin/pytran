# pyhitran

A python package to read and use the HITRAN database.

### Intro
This is my attempt to create a general interface into the HITRAN
database for the purpose of plotting and using IR absorption cross
sections. I also use the files provided in the 'Global Data' folder to
aid in correctly evaluating the cross sections.

### Design Goals

##### hitran.py
- [x] Read HITRAN .par files into line list dictionary
- [x] Compute cross section from line list
- [x ] Use partition function in computing cross section
- [ ] Create utilities to interpret line lists (i.e. list of
  transitions with bounds and strengths, compute band Einstein
  coefficient, etc)
- [ ] Determine HITRAN_ROOT from shell (eg where are the .par files
located)
- [ ] make 'get_filename' that take mol number or name and give HITRAN
.par file

##### hitran_supplemental.py
- [x] Read HITRAN molparam.txt file into molecule_param dictionary
- [x] Create utilities to return names, indicies, values for molparam
dict
- [x] Read HITRAN parsum.dat into parsum dictionary
- [x] Create utility to evaluate parsum for particular isotopologue
- [ ] Implement low temperature partition sum

##### Long term goal
- [ ] Develop into a Python Package
- [ ] Make a GUI like pyHAWKS

### Aknowledgements

I adapted the initial code from https://github.com/nzhagen/hitran,
thanks go to the author.

