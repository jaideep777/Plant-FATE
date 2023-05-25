#include <vector>
#include <ostream>

#include "plant_params.h"
#include "plant_geometry.h"
#include "assimilation.h"
#include "plant.h"

#include "climate.h"
#include "light_environment.h"


class ErgodicEnvironment : public env::Climate, public LightEnvironment {
	public:
	void print(double t);
};


class LifeHistoryOptimizer{
	public:
	plant::Plant P;
	ErgodicEnvironment C;

	std::string params_file;
	std::string met_file = "";
	std::string co2_file = "";

	double dt = 0.1; 

	double rep;
	double litter_pool;
	double seeds;
	double prod;

	public:

	void init();

	void printHeader(std::ostream &lfout);

	void printState(double t, std::ostream& lfout);
	void printPlant();

	void set_traits(std::vector<double> tvec);
	std::vector<double> get_traits();

	void set_state(std::vector<double>::iterator it);

	void get_rates(std::vector<double>::iterator it);

	void grow_for_dt(double t, double dt);


	double calcFitness();

};


