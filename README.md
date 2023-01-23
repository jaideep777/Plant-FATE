# PlantFATE Eco-evolutionary vegetation model

Plant-FATE, standing for Plant FunctionAl Trait Evolution, is an eco-evolutionary vegetation model that accounts for multi-timescale adaptations of invdividual plants and plant species to the environment.

## Prerequisites

Plant-FATE can be built without any third party libraries as shown here (but GSL and Eigen can be optionally used - we will not cover these cases for now). 

You do need git to select the appropriate branches for compiling. Instructions to install git can be found here: https://github.com/git-guides/install-git

## Installation

You need three repos to run PlantFATE: 
1) [P-hydro](https://github.com/jaideep777/phydro): photosynthesis model
2) [libpspm](https://github.com/jaideep777/libpspm): PSPM solver library
3) [PlantFATE](https://github.com/jaideep777/Plant-FATE): The PlantFATE EVM

Create a folder that will serve as a root folder for cloning these repos. Clone all three repos in same (root) folder, so that your directory structure looks like this:

```
  root
  |---- phydro
  |     |--- inst/include
  |
  |---- libpspm
  |     |--- include
  |     |--- lib
  |
  |---- Plant-FATE
  |     |--- inst/include

```

For each repo, we clone it and checkout the appropriate branch. The latest branches that should be used are shown in the commands below. 

```
cd <path/to/root>

git clone https://github.com/jaideep777/phydro
cd phydro
git checkout master
cd ..

git clone https://github.com/jaideep777/libpspm
cd libpspm
git checkout develop
cd ..

git clone https://github.com/jaideep777/Plant-FATE
cd Plant-FATE
git checkout feature_evolution
cd ..

```

### Prepare input data and simulation configuration

1. Paste data files in folder tests/data

2. Simulation config is in `tests/params/p.ini`. Change these parameters as necessary. Particularly, set the output folder: output files will be generated in the folder specified in a folder hierarchy as follows:  Plant_FATE/`outDir`/`exptName`


### R Installation

#### Prerequisites 

If you are using R in Windows, you also need to install a C++ compiler for R and the Rcpp package. Please follow these guidelines to install both: https://teuder.github.io/rcpp4everyone_en/020_install.html  

#### Compiling the PlantFATE R package

After cloning the three repos and checking out the appropriate branches, create a new R Project, select 'Exisiting Directory' in the wizard, and select the Plant-FATE folder. At the end of this step, you should have the following files in the Plant-FATE folder (among others):

```
  root
  |---- Plant-FATE
  |     |--- ...
  |     |--- Plant-FATE.Rproj
  |     |--- DESCRIPTION
  |     |--- NAMESPACE
  |     |--- vignettes
  |     |--- ...
```

Open the project and simply perform a build (type ctrl+shift+B or via the menu, click on Build/Install Package). In Windows, if you get a prompt saying you need additional build tools, go ahead and install them.

Please see the vignettes folder for demos on how to use the PlantFATE and TreeLife models in R.

### Native C++ Installation

Perform these steps from the root folder:

#### Compile and test `libpspm`

```
cd libpspm
make clean testclean
make
make check
```

*(all tests should pass)*

#### Test `Phydro` 

```
cd phydro
make clean testclean
make check
```

*(all tests should pass)*

#### Compile and run PlantFATE

The main program to run Plant-FATE is in `tests/pf.cpp`:

```
make clean testclean
make check TEST_FILES=tests/pf.cpp
```

## Analyse results

Analyse results with `tests/Rscripts/amazon_calibration.R` 
Dont forget to change path to Working dir

## Author Contact

Jaideep Joshi (jaideep777@gmail.com)

## Acknowledgements

This project has received funding from the European Union’s Horizon 2020 research and innovation programme under the Marie Skłodowska-Curie Actions fellowship (grant agreement No 841283).

