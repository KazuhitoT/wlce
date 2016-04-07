# WLCE code

Wang-Landau sampling with Cluster Expansion

## Description
version 0.6

***DEMO for two-dimentional square lattice Ising spin:***

```
cd 2dising
../bin/cluster -d=1.0
cp multiplicity.out multiplicity.in
cp clusters.out clusters.in
cp labels.out labels.in
../bin/wang-landau
```

## Features
- Canonical Metropolis sampling
- Canonical Wang-Landau sampling
- Semi-Grand-Canonical  Wang-Landau sampling

## Requirement

- C++11 compiler
- CMake
- Boost
- Eigen

Confirmed
- gcc >= version 4.8.4
- clang version 7.0.2
- icpc version 16.0.1

## Installation

```
git clone https://github.com/kazuhitoT/wlce
cmake .
make
```

## Usage

### cluster
input file
- poscar.in

output file
- clusters.out, multiplicity.out, labels.out

command line option
- -d=[truncation distance]
- -p=[filename for poscar.in]

### corrdump
input file
- corrdump.ini, clusters.in, multiplicity.in, labels.in

output file
- none

corrdump.ini
```
[INPUT]
SPINCE = -1 1
```

### metropolis
input file
- metropolis.ini, clusters.in, multiplicity.in, labels.in, ecicar

output file
- metropolis.out

corrdump.ini
```
[INPUT]
SPINCE        = -1 1
SPINPOSCAR    = -1 1
MCSTEP        = 10
SAMPLESTEP    = 10
TEMPERATURE   = 700
or TEMPERATURE   = 900 100 100
or TEMPERATURE   = 100 100 900
```

### wang-landau
input file
-   wang-landau.ini, poscar.spin, ecicar, poscar.in(same with cluster), clusters.in(output by cluster), multiplicity.in(output by cluster), labels.in, ecicar

output file
- out-wl*.dat, macrostate.out

wang-landau.ini
```
[INPUT]
SPINCE        = -1 1
SPINPOSCAR    = -1 1
FLATCHECKSTEP = 100
MCSTEP        = 100
LOGFACTOR     = 4
LOGFLIMIT     = 0.0000001
BIN           = 255
EMIN          = -1
EMAX          = 1
FLATCRITERION = 0.8
CHEMIPOT      = 0 1 # sgc mode
```


## Author

Kazuhito TAKEUCHI
[@kazuhitoT](https://github.com/kazuhitoT)

## License
[MIT](https://opensource.org/licenses/mit-license.php)


## For more details on Cluster Expansion
> * J.M. Sanchez et. al., Physica A 128, 334-350 (1984).
> * A. van de Walle and M. Asta, Modelling Simul. Mater. Sci. Eng. 10, 521-538 (2002).


## For more details on Wang-Landau sampling
> * F. G. Wang and D. P. Landau, Phys. Rev. Lett. 86, 2050 (2001).
> * F. G. Wang and D. P. Landau, Phys. Rev. E 64, 056101 (2001).
