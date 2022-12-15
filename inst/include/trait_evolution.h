#ifndef PLANT_FATE_TRAIT_EVOLUTION_H_
#define PLANT_FATE_TRAIT_EVOLUTION_H_

#include <solver.h>
#include <vector>
#include <string>
#include "utils/moving_average.h"

// Extend the Species class from libpspm to allow trait evolution
template <class Model>
class MySpecies : public Species<Model>{
	public:
	std::string species_name;
	bool isResident = true;
	double t_introduction = 0;

	double fg_dx = 0.001;
	double invasion_fitness;
	double r0;
	std::vector<MySpecies<Model>*> probes;
	std::vector<double> fitness_gradient;
	std::vector<double> trait_variance;
	std::vector<double> trait_scalars;     // these scalars will be applied to fg_dx
	std::vector<std::string> trait_names;
	
	MovingAverager seeds_hist;
	MovingAverager r0_hist;

	MySpecies(Model M, bool res=true);

	void set_traits(std::vector<double> tvec);
	std::vector<double> get_traits();

	// Species(vector<double> tvec, double u0, bool res);
	
	void createVariants(Model M);
	
	void calcFitnessGradient();
	void evolveTraits(double dt);

	void print_extra();
};

#include "trait_evolution.tpp"

#endif

