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
	ErgodicEnvironment();
	void print(double t);
};


class LifeHistoryOptimizer{
	public:
	plant::Plant P;
	ErgodicEnvironment C;

	// std::string paramsFile;
	// std::string met_file = "";
	// std::string co2_file = "";

	plant::PlantParameters par0;
	plant::PlantTraits traits0;

	io::Initializer I;

	double dt = 0.1; 

	double rep;
	double litter_pool;
	double seeds;
	double prod;

	public:

	LifeHistoryOptimizer(std::string params_file);
	
	void set_metFile(std::string metfile);
	void set_co2File(std::string co2file);
	
	void init();

	std::vector<std::string> get_header();
	void printHeader(std::ostream &lfout);

	std::vector<double> get_state(double t);
	void printState(double t, std::ostream& lfout);

	void printMeta();

	void set_traits(std::vector<double> tvec);
	std::vector<double> get_traits();

	void set_state(std::vector<double>::iterator it);

	void get_rates(std::vector<double>::iterator it);

	void grow_for_dt(double t, double dt);


	double calcFitness();

};


