#ifndef PLANT_FATE_ENV_LIGHT_ENVIRONMENT_H_
#define PLANT_FATE_ENV_LIGHT_ENVIRONMENT_H_

#include <iostream>
#include <vector>

namespace env{

class LightEnvironment {
	public:
	bool use_ppa = false;

	// Smooth environment 

	// PPA environment	
	int n_layers;
	std::vector<double> z_star;
	std::vector<double> canopy_openness;
	
	
	public:
	LightEnvironment();

	//double LightEnvironment::LA_above_z(double t, double z, Solver *S);

	//// This function must do any necessary precomputations to facilitate evalEnv()
	//// Therefore, this should calculate env for all X when it is a function of X
	//// In such a case, the solver's SubdivisionSpline can be ussed
	//// Note: The state vector in the solver will not be updated until the RK step is completed. 
	//// Hence, explicitly pass the state to this function.
	//// ~
	//// Also this is the only function that exposes the state vector, so if desired, the state vector 
	//// can be saved from here and reused in other rate functions (using createIterators_state())
	//// ~
	//// TODO: In Solver, add a add_iAttribute() function, that will calculate some individual 
	//// level attributes from x, which can be reused if required. E.g., in Plant, we can add leaf_area
	//// as an iAttribute. iAttributes can be mapped to integers, say using enums
	//// Alternatively, switch to Indiviudual class as a template parameter for solver
	//void LightEnvironment::computeEnv(double t, Solver * S);

	void print();
	
};


} // env

#endif
