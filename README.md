# PlantFATE Eco-evolutionary vegetation model

Plant-FATE, standing for Plant FunctionAl Trait Evolution, is an eco-evolutionary vegetation model that accounts for multi-timescale adaptations of invdividual plants and plant species to the environment.

## Prerequisites

### Prepare input data and simulation configuration

1. Simulation configuration should be specified in `.ini` format. A template can be found in `tests/params/p_test_v2.ini`. Change these parameters as necessary. 

2. Paste data files in folder tests/data. The format for input data can be found in the comments in the template ini file.


## Installation

### Native C++ Installation

#### Compilation

To compile PlantFATE, navigate to the downloaded directory and issue the following commands in a linux terminal: 

```
make
make check
```

#### Model run

The main executables to run PlantFATE and TreeLife are in `Plant-FATE/bin`. It may be convenient to add this path to your `~/.bashrc` file. To do so, add the following line to your bashrc:
```
export PATH=$PATH:<path/to/Plant-FATE>/bin
```

To run Plant-FATE, use the following syntax:
```
plantfate <path/to/p.ini> start_year end_year
```

To run TreeLife, use the following syntax:
```
treelife <path/to/p.ini> n_years
```


### R Installation

#### Prerequisites 

If you are using R in Windows, you also need to install a C++ compiler for R and the Rcpp package. Please follow these guidelines to install both: https://teuder.github.io/rcpp4everyone_en/020_install.html  

You also need to install the devtools package to install Plant-FATE directly from Github. This can be installed using the following command:
```
install.packages("devtools")
```

#### Installing the PlantFATE R package

You can install Plant-FATE directly from github with the following command:

```
devtools::install_github("jaideep777/Plant-FATE", ref="develop", force = T)
```

Please see the vignettes folder for demos on how to use the PlantFATE and TreeLife models in R.


### Python installation

Python3 and pip3 is required to install PlantFATE as a python library. Run: 

```
make python
```


## Analyse results

Analyse results with `tests/Rscripts/pf_test_analysis.R` 

Dont forget to change path to Working dir

## Lead author Contact

Jaideep Joshi (jaideep777@gmail.com)

For all contributors, see `AUTHORS.md`

## Acknowledgements

This project has received funding from the European Union’s Horizon 2020 research and innovation programme under the Marie Skłodowska-Curie Actions fellowship (grant agreement No 841283).

