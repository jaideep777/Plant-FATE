#ifndef PLANT_FATE_TRAIT_EVOLUTION_H_
#define PLANT_FATE_TRAIT_EVOLUTION_H_

#include <solver.h>
#include <vector>

// Extend the Species class from libpspm to allow trait evolution
template <class Model>
class MySpecies : public Species<Model>{
	public:
	bool isResident = true;

	double fg_dx = 0.001;
	double invasion_fitness;
	double r0;
	std::vector<Species<Model>*> probes;
	std::vector<double> fitness_gradient;
	std::vector<double> trait_variance;
	std::vector<double> trait_scalars;     // these scalars will be applied to fg_dx

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

