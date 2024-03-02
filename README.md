# PlantFATE Eco-evolutionary vegetation model

Plant-FATE, which stands for Plant Functional Acclimation and Trait Evolution, is an eco-evolutionary vegetation model designed to predict the multi-timescale adaptations of individual plants, plant species, and mixed-species communities to changing environmental conditions.

## Prerequisites

### Prepare Input Data and Simulation Configuration

1. **Simulation Configuration**: The simulation configuration should be specified in `.ini` format. You can find a template in `tests/params/p_test_v2.ini`. This file contains parameters that control various aspects of the simulation such as plant traits, environmental conditions, and simulation duration. Modify these parameters as necessary according to your simulation requirements.

2. **Input Data**: Specify the paths to your data files in your configuration file. The format for input data can be understood from the comments within the template configuration file. 

## Installation

### Native C++ Installation

Native C++ installation is currently supported only on Linux machines.

#### Compilation

To compile PlantFATE from source code, follow these steps:

1. Navigate to the downloaded directory containing the PlantFATE source code.
2. Open a terminal in Linux.
3. Run the following commands:
   ```bash
   make
   make check
   ```
   These commands will compile the PlantFATE source code and run tests to ensure that the compilation was successful.

#### Model Run

1. The main executables to run PlantFATE and TreeLife are located in the `bin` directory. You may want to add this path to your `~/.bashrc` file for convenience. To do so, add the following line to your bashrc:
   ```bash
   export PATH=$PATH:<path/to/Plant-FATE>/bin
   ```
   Replace `<path/to/Plant-FATE>` with the actual path to the `Plant-FATE` directory on your system.

2. To run PlantFATE, use the following syntax:
   ```bash
   plantfate <path/to/p.ini> start_year end_year
   ```
   Replace `<path/to/p.ini>` with the path to your simulation configuration file (e.g., `p_test_v2.ini`), `start_year` with the starting year of the simulation, and `end_year` with the ending year.

3. To run TreeLife, use the following syntax:
   ```bash
   treelife <path/to/p.ini> n_years
   ```
   Replace `<path/to/p.ini>` with the path to your simulation configuration file (e.g., `p_test_v2.ini`), and `n_years` with the number of years over which lifetime fitness is calculated. The argument `n_years` is optional and defaults to 500 years.

### R Installation

#### Prerequisites

1. **C++ Compiler**: If you're using R on Windows, you need to install a C++17-compatible compiler for R and the Rcpp package. Follow these guidelines to install both: [Rcpp Installation Guide](https://teuder.github.io/rcpp4everyone_en/020_install.html).

2. **Devtools Package**: Install the devtools package using the following command in R:
   ```R
   install.packages("devtools")
   ```

#### Installing the PlantFATE R Package

You can install Plant-FATE directly from GitHub using the devtools package in R. Run the following command:
```R
devtools::install_github("jaideep777/Plant-FATE", ref="develop", force = T)
```

Please see the vignettes folder for demos on how to use the PlantFATE and TreeLife models in R.

### Python Installation

Python3 and pip3 are required to install PlantFATE as a Python library. Follow these steps:

1. Open a terminal.
2. Run the following command:
   ```bash
   make python
   ```
   This command will install PlantFATE as a Python library.

## Analyze Results

You can analyze the results using the `tests/Rscripts/pf_test_analysis.R` script. Ensure to change the path to the working directory as necessary.

## Lead Author Contact

For inquiries, you can contact the lead author, Jaideep Joshi, at jaideep777@gmail.com.

For a full list of contributors, please refer to `AUTHORS.md`.

## Acknowledgments

This project has received funding from the European Union’s Horizon 2020 research and innovation program under the Marie Skłodowska-Curie Actions fellowship (grant agreement No. `841283`) and from the Strategic Initiatives Program of IIASA (project `RESIST`).

