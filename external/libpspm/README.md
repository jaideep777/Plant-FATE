---
output:
  word_document: default
  html_document: default
---
# libpspm

libpspm is a C++ library to simulate Physiologically Structured Population Models (PSPMs).

It implements four numerical methods to discretize the McKendrick-von Foerster Partial Differential Equation (PDE), which forms the core of a PSPM, into a set of Ordinary Differential Equations (ODEs); see Zhang et al. 2017:

1) Fixed Mesh Upwind (FMU)
2) Implicit Fixed Mesh Upwind (IFMU) 
3) Escalator Boxcar Train (EBT)
4) Characteristic Method (CM)

The resulting ODEs can be solved using different ODE-solvers. Two such are currently implemented:

1) Runge-Kutta 4-5 Cash-Karp
2) LSODA

libpspm supports multi-species systems with multiple size-structured and unstructured species.


## Installation

To compile as a static library, just perform 
```
make 
make check
```
This will create a static library (`libpspm.a`) in the `lib` folder. Link to this library as `-lpspm` while compiling your program.

Alternatively, just include the contents of the `src` and `include` folders in your respective project folders, compile each `cpp` file into an object file, and link with the objects of your project. 

### Stable release

Coming soon.

### Development version

Latest development version can be found here: https://github.com/jaideep777/pspm/tree/develop


## Usage

To solve a PSPM, you need to define two classes for the `Environment` and the `Individual` (not necessarily with those names). The `Environment` class must be derived from the `EnvironmentBase` class provided by the library. Although the `Individual` class can be defined as a standalone class, the easiest way to define it is also to  derive it from the `IndividualBase` class provided by the library. In these derived classes you must implement the following member functions (one in `Environment` and four in `Individual`):

```C++
class Environment : public EnvironmentBase{
	public:
	// This function should perform environment computation, given the current state 
	void computeEnv(double t, Solver * S, std::vector<double>::iterator s, std::vector<double>::iterator dsdt){
	}
};

class Individual : public IndividualBase{
	public:
	// return the initial density (initial condition) as a function of 
	// the physiological variable x, environment _env, and input flux of newborns bf
	double init_density(double x, void * _env, double bf){
	}

	// return the growth rate as a function of 
	// the physiological variable x, time t, and environment _env
	double growthRate(double x, double t, void * _env){
	}

	// return the mortality rate as a function of 
	// the physiological variable x, time t, and environment _env
	double mortalityRate(double x, double t, void * _env){
	}

	// return the fecundity (birth rate) as a function of 
	// the physiological variable x, time t, and environment _env
	double birthRate(double x, double t, void * _env){
	}

};
```

Having defined these classes, create an `Environment` object and a species of `Individual`s as follows:

```C++
Environment E;
Species<Individual> spp;
```

Then create a `Solver` by specifying the PSPM method and the ODE solver, e.g.,

```C++
Solver S(SOLVER_FMU, "rk45ck");
```

Add the species and the environment to the solver:
```C++
// add the created Species 'spp' to the solver along with some species-specific properties 
S.addSpecies(25, xb, xm, false, &spp, 0, -1);
S.setEnvironment(&E);
```

Initialize and run the PSPM:
```C++
S.resetState();
S.initialize();
for (double t=0; t < 10; t=t+1) {
    S.step_to(t);
}
```

Concrete examples can be found in the `demos` folder. For more detailed tutorials explaining the demos, see:

1. Daphnia model - 
2. RED model - 

## Author and contact

Jaideep Joshi
joshi@iiasa.ac.at


## Acknowledgement

This project was funded by a Marie Sklodowska-Curie fellowship H2020-MSCA-IF-2019, project PlantFATE.



