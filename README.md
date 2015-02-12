# pyhitran

A python package to read and use the HITRAN database.

### Intro
This is my attempt to create a general interface into the HITRAN
database for the purpose of plotting and using IR absorption cross
sections. I also use the files provided in the 'Global Data' folder to
aid in correctly evaluating the cross sections.

### Design Goals
- [x] Read HITRAN .par files into line list dictionary
- [x] Compute cross section from line list
- [ ] Use partition function in computing cross section
- [ ] Create utilities to interpret line lists (i.e. list of
  transitions with bounds and strengths, compute band Einstein
  coefficient, etc)
- [x] Read HITRAN molparam.txt file into molecule_param dictionary
- [ ] Create utilities to return names, indicies, values for molparam
dict
- [x] Read HITRAN parsum.dat into parsum dictionary
- [] Create utility to evaluate parsum for particular isotopologue

### Aknowledgements



