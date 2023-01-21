#include <vector>
#include <ostream>

#include "plant_params.h"
#include "plant_geometry.h"
#include "assimilation.h"
#include "plant.h"

#include "climate.h"
#include "light_environment.h"


class ErgodicEnvironment : public env::Climate, public env::LightEnvironment {
	public:
	void print(double t);
};


class LifeHistoryOptimizer{
	public:
	plant::Plant P;
	ErgodicEnvironment C;

	double dt = 0.1; 

	double rep;
	double litter_pool;
	double seeds;
	double prod;

	public:

	void init();

	void printHeader(std::ostream &lfout);

	void printState(double t, std::ostream& lfout);

	void set_state(std::vector<double>::iterator it);

	void get_rates(std::vector<double>::iterator it);

	void grow_for_dt(double t, double dt);


	double calcFitness();

};


