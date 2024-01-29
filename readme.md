## Autocorrelation Analysis Tools Accelerated for LAMMPS Trajectory


## install

```shell
# clone the source codes
git clone https://github.com/zhijs/autocorrelation-for-lammps.git

```

## compile

### fortran code
```shell
# compile the fortran lib
source xxxx/setvars.sh 

ifx -O3 -qopenmp -o cal_neighbour.f90 ./cal_neighbour.f90.f90 

ifx -O3 -qopenmp -o cal_neighbour_for_LLPS.f90 ./cal_neighbour_for_LLPS.f90
```

### C++ code
```shell
# compile the cpp code
c++ inside_ac.cpp -std=c++17 -O3 -o inside_ac
c++ exchange.cpp -std=c++17 -O3 -o exchange
c++ full_ac.cpp -std=c++17 -O3 -fopenmp -o  full_ac # only full_ac supports openmp to accelerate
```

## Usage

### First step: produce neighbour list files

for example: if you want to calculate the autocorrelation that TF binding chromatin, run followed command will product neighbour.txt file

```shell
# bead.dat and polymer.dat represent the TF data and chromatin data repectively
./cal_neighbour bead.dat polymer.dat
```


### Second step: run autocorrelation analysis

1. Calculate autocorrelation between different atom types through a LAMMPS trajectory

```shell
./full_ac -i './data/polymer_1.dat' -id '1_2' -n './data/neighbors.txt' -t 1000  -c 5 -freq 500 -si 0 -ei 8000000 -o './data/full_ac.txt'
```

2. Calculate autocorrelation between different atom types through a LAMMPS trajectory when phase-phase separation occurs

```shell
# firstly, calculate condenstate cluster and the exchanged rate of condenstate can also be calculated
 ./exchange -i './data/polymer_1.dat' -cf './data0.65/cluster.txt' -n './origin/neighbors.txt' -t 1000  -c 5 -freq 500 -bw 80 -bh 80 -bl 80 -o './data/output.txt'

# Secondly, run the autocorrelation program and make use of the cluster result
./inside_ac -cf './data0.53/cluster/cluster_1.txt' -n "./data/neighbors.txt" -id '1_2' -t 1000 -freq 500 -o './data/cluster_ac.txt' 

```

Parameter Description:

* i: 
LAMMPS data file, it can be the first data file
* id: 
AtomA_AtomB, use '_' symbol to split the different atom types which you want to analyze autocorrelation
* n:
LAMMPS neighbor list file
* t:
Frame number
* c:
CutOff, for the LLPS cluster computation
* freq: 
Timestep of per frame
* -si: 
(in order to save computing time) Atom start index, recommend setting 0
* -ei: 
(in order to save computing time) Atom end index, recommend setting maximum atom index


