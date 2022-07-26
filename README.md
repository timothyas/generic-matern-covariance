# A Practical Formulation for an Anisotropic and Nonstationary Matérn Class Correlation Operator

[![DOI](https://zenodo.org/badge/471468580.svg)](https://zenodo.org/badge/latestdoi/471468580)

## TODO

- [ ] MITgcm code to separate branch
- [ ] Clean this repo
    - move relevant pych scripts here
    - clean up / organize
    - rosy pig ... public and rename
- [ ] Document rotation case

## MITgcm TODO


- [ ] First step...
    - [x] Reorganize WC01 code
    - [x] Make 2D and 3D num operators separate...
    - [x] Change maskname from 'maskC' to just 'C'
    - [ ] Test compile / setup
    - [ ] update smooth check and printing in readparms

- [ ] Initial port of 2D/3D routines
    - [x] initial import
    - [x] Incorporate 3D updates (e.g. tolerance, num applications) to 2D fields
    - [x] Rename jacobi -> SOR
    - [ ] Test compile / setup
    - [ ] update smooth check and printing in readparms
    - [ ] add redi stuff?
    - [ ] Documentation

- [ ] 2D open boundary XZ / YZ routines
    - smoothdims character arg
- [ ] eigenvector decomposition routines
