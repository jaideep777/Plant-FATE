#include "pspm_environment.h"
#include "pspm_plant.h"
#include <species.h>



LightEnvironment::LightEnvironment(double openness) : Environment(openness){
}

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
void LightEnvironment::computeEnv(double t, Solver * S, std::vector<double>::iterator s, std::vector<double>::iterator dsdt){
	//            _xm 
	// Calculate / w(z,t)u(z,t)dz
	//        xb`
	auto canopy_openness = [S, t, this](double z){
		double kI = 0.5;

		double leaf_area_above_z = 0;
		
		// Loop over resident species --->
		for (int k=0; k<S->species_vec.size(); ++k){
			auto la_above = [z,k,S](int i, double t){
				auto& p = ((Species<PSPM_Plant>*)S->species_vec[k])->getCohort(i);
				double a = p.area_leaf_above(z, p.vars.height, p.vars.area_leaf);
//				std::cout << "(" << i << "," << a << ")" << "\t";
				return a;	
			};
			leaf_area_above_z += S->integrate_wudx_above(la_above, t, {z}, k);
			//leaf_area_above_z += S->integrate_x(la_above, t, state_vec, i);
		}

//		std::cout << "LA(" << z << ") = " << exp(-kI*leaf_area_above_z) << "\n";
		return exp(-kI*leaf_area_above_z);
	};	

	//cout << S->xb << " " << S->getMaxSize() << endl;	
	//time = t;
	//for (int s=0; s<S->n_species(); ++s) S->get_species(s)->u0_save = S->get_u0(t, s);
	double max_size = 0;
	for (int k=0; k<S->species_vec.size(); ++k){
		max_size = std::max(max_size, S->maxState(k)[0]);
	}
	light_profile.construct(canopy_openness, 0, max_size);
}






FixedEnvironment::FixedEnvironment(double o){
	openness = o;
}

double FixedEnvironment::canopy_openness(double z) const {
	return openness;
}

void FixedEnvironment::computeEnv(double t, Solver * S, std::vector<double>::iterator s, std::vector<double>::iterator dsdt){
}

double FixedEnvironment::patch_survival(double t) const{
	return 1;
}

