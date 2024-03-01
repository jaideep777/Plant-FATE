#ifndef PLANT_FATE_PFATE_LIFEHISTORY_OPTIMIZER_H_
#define PLANT_FATE_PFATE_LIFEHISTORY_OPTIMIZER_H_

#include <vector>
#include <ostream>

#include "traits_params.h"
#include "plant_architecture.h"
#include "assimilation.h"
#include "plant.h"

#include "climate.h"
#include "climate_stream.h"
#include "light_environment.h"

namespace pfate{

class ErgodicEnvironment : public env::LightEnvironment, public env::Climate {
	public:
	ErgodicEnvironment();
	void print(double t) override;

	// override computeEnv() to NOT update the light profile
	void computeEnv(double t, Solver * sol, std::vector<double>::iterator S, std::vector<double>::iterator dSdt) override;
};


class LifeHistoryOptimizer{
	public:
	plant::Plant P;
	ErgodicEnvironment C;
	env::ClimateStream c_stream;

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
	
	void set_i_metFile(std::string file);
	void set_a_metFile(std::string file);
	void set_co2File(std::string co2file);
	void init_co2(double _co2);

	void init();

	void update_climate(double julian_time);

	std::vector<std::string> get_header();
	void printHeader(std::ostream &lfout);

	std::vector<double> get_state(double t);
	void printState(double t, std::ostream& lfout);

	void printMeta();

	// void set_traits(std::vector<double> tvec);
	// std::vector<double> get_traits();

	void set_state(std::vector<double>::iterator it);

	void get_rates(std::vector<double>::iterator it);

	void grow_for_dt(double t, double dt);


	double calcFitness();

};

} // namespace pfate

#endif
