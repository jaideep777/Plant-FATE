#ifndef PLANT_FATE_TRAIT_EVOLUTION_H_
#define PLANT_FATE_TRAIT_EVOLUTION_H_

#include <solver.h>

// Extend the Species class from libpspm to allow trait evolution
template <class Model>
class MySpecies : public Species<Model>{
	public:
	double r0;

	MySpecies(Model M);
};

#include "trait_evolution.tpp"

#endif

