# WLCE code

Wang-Landau sampling with Cluster Expansion

## Description
version 0.7

***DEMO for two-dimentional square lattice Ising spin:***

```
cd 2dising
../bin/cluster -d=1.0
cp multiplicity.out multiplicity.in
cp clusters.out clusters.in
cp labels.out labels.in
../bin/wang-landau

gnuplot
>plot "out-wl20.dat" u ((($2)+($3))/2.):4
```

format for output file (out-wl*.dat)
```
# index bin_energy_min bin_energy_max entropy histogram
0 -1.0000000000 -0.9921568627 30340.0000000000 7585.0000000000
1 -0.9921568627 -0.9843137255 30412.0000000000 7603.0000000000
2 -0.9843137255 -0.9764705882 30380.0000000000 7595.0000000000
3 -0.9764705882 -0.9686274510 30784.0000000000 7696.0000000000
4 -0.9686274510 -0.9607843137 30884.0000000000 7721.0000000000
5 -0.9607843137 -0.9529411765 30916.0000000000 7729.0000000000
6 -0.9529411765 -0.9450980392 31012.0000000000 7753.0000000000
7 -0.9450980392 -0.9372549020 31100.0000000000 7775.0000000000
8 -0.9372549020 -0.9294117647 31252.0000000000 7813.0000000000
9 -0.9294117647 -0.9215686275 31272.0000000000 7818.0000000000
...
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
cd wlce
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
- -p=[filename (default: poscar.in)]

### getconf
input file
- getconf.ini, poscar.in clusters.in, multiplicity.in, labels.in,

output file
- none

getconf.ini
```
[INPUT]
SPINCE     = -1 1
SPINPOSCAR = -1 1
```

command line option
- -d=[truncation distance] (same in cluster)
- none output total energy and compositions (need ecicar)
- -c output correlation functions

### metropolis
input file
- poscar.spin, metropolis.ini, clusters.in, multiplicity.in, labels.in, ecicar

output file
- metropolis.out

metropolis.ini
```
[INPUT]
SPINCE        = -1 1
SPINPOSCAR    = -1 1
MCSTEP        = 10
SAMPLESTEP    = 10
TEMPERATURE   = 700
or TEMPERATURE   = 900 100 100 # temperature decreases from 900 to 100 with 100 step
or TEMPERATURE   = 100 100 900 # temperature increases vice versa.
```

format for output file (metropolis.out)
```
# temperature internal_energy variance_of_internal_energy
900.000000 -0.012319 0.000688
800.000000 -0.021950 0.000480
700.000000 -0.006619 0.000491
600.000000 -0.016833 0.000556
500.000000 -0.012944 0.000580
400.000000 -0.017795 0.000536
300.000000 -0.019312 0.000498
200.000000 -0.035830 0.000739
...
```

### wang-landau
input file
-   wang-landau.ini, poscar.spin, ecicar, clusters.in(output by cluster), multiplicity.in(output by cluster), labels.in, ecicar

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
