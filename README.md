# PlantFATE Eco-evolutionary vegetation model

Plant-FATE, standing for Plant FunctionAl Trait Evolution, is an eco-evolutionary vegetation model that accounts for multi-timescale adaptations of invdividual plants and plant species to the environment.

## Prerequisites

You need git, gsl, and eigen to be able to run Plant-FATE

```
sudo apt-get install git libeigen3-dev libgsl-dev
```

## Installation

You need three repos to run PlantFATE: 
1) [P-hydro](https://github.com/jaideep777/phydro): photosynthesis model
2) [libpspm](https://github.com/jaideep777/libpspm): PSPM solver library
3) [PlantFATE](https://github.com/jaideep777/Plant-FATE): The PlantFATE EVM

You can clone all of them in the same root folder

### Install and test `libpspm`

```
git clone https://github.com/jaideep777/libpspm
git checkout hotfix_afterstep
make clean testclean
make
make check
```

**(all tests should pass)**

### Install and test `Phydro` 

```
git clone https://github.com/jaideep777/phydro
git checkout master
make clean testclean
make check
```

**(all tests should pass)**

### Set up PlantFATE

To set up PlantFATE, follow these steps:

1. clone the repo and checkout the desired branch.

```
git clone https://github.com/jaideep777/Plant-FATE
git checkout feature_evolution
```

2. Change paths in the makefile to correctly point to eigen, phydro, and libpspm

3. Paste data files in folder tests/data

4. Change parameters in tests/params/p.ini (if necessary)

5. Compile and run:
```
make clean testclean
make check TEST_FILES=tests/plantfate_evolution_test.cpp
```
### Analyse results

Analyse results with `tests/Rscripts/amazon_calibration.R` 
Dont forget to change path to Working dir

## Author Contact

Jaideep Joshi (jaideep777@gmail.com)

## Acknowledgements

This project has received funding from the European Union’s Horizon 2020 research and innovation programme under the Marie Skłodowska-Curie Actions fellowship (grant agreement No 841283).

