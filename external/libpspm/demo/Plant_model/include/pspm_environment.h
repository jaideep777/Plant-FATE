#ifndef PLANT_PSPM_ENV_H_
#define PLANT_PSPM_ENV_H_

#include "environment.h"
#include <environment_base.h>
#include <solver.h>

class LightEnvironment : public plant::Environment, public EnvironmentBase {
	//double evalEnv(double x, double t){
	//    env.light_profile.eval(x); // return 1;
	//}

	public:
	LightEnvironment(double openness);

	// This function must do any necessary precomputations to facilitate evalEnv()
	// Therefore, this should calculate env for all X when it is a function of X
	// In such a case, the solver's SubdivisionSpline can be ussed
	// Note: The state vector in the solver will not be updated until the RK step is completed. 
	// Hence, explicitly pass the state to this function.
	// ~
	// Also this is the only function that exposes the state vector, so if desired, the state vector 
	// can be saved from here and reused in other rate functions (using createIterators_state())
	// ~
	// TODO: In Solver, add a add_iAttribute() function, that will calculate some individual 
	// level attributes from x, which can be reused if required. E.g., in Plant, we can add leaf_area
	// as an iAttribute. iAttributes can be mapped to integers, say using enums
	// Alternatively, switch to Indiviudual class as a template parameter for solver
	void computeEnv(double t, Solver * S, std::vector<double>::iterator s, std::vector<double>::iterator dsdt);


};


class FixedEnvironment : public EnvironmentBase {
	private:
	double openness;

	public:
	FixedEnvironment(double o);

	double canopy_openness(double z) const;

	void computeEnv(double t, Solver * S, std::vector<double>::iterator s, std::vector<double>::iterator dsdt);

	double patch_survival(double t) const;

	
};

#endif
