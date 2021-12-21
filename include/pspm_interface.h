#ifndef PLANT_FATE_PSPM_INTERFACE_H_
#define PLANT_FATE_PSPM_INTERFACE_H_

#include <solver.h>
#include "light_environment.h"
#include "climate.h"
#include "plant.h"


class PSPM_Plant : public plant::Plant {
	public:
	
	double t_birth = 0;

	std::vector<std::string> varnames = {"name", "|lma|", "D", "g", "lai", "mort", "seeds"}; // header corresponding to the print function below
	std::vector<std::string> statevarnames = {"lai", "mort", "seedpool"};                 // header corresponding to state output

	int nrc = 0; // number of evals of compute_vars_phys() - derivative computations actually done by plant
	int ndc = 0; // number of evals of mortality_rate() - derivative computations requested by solver
	int nbc = 0;

	PSPM_Plant(); 

	void set_size(double _x);

	double init_density(double x, void * _env, double input_seed_rain);

	void preCompute(double x, double t, void * _env);

	double establishmentProbability(double t, void * _env);

	double growthRate(double x, double t, void * _env);
	double mortalityRate(double x, double t, void * _env);
	double birthRate(double x, double t, void * _env);
	
	void init_state(double t, void * _env);

	std::vector<double>::iterator set_state(std::vector<double>::iterator &it);

	std::vector<double>::iterator get_state(std::vector<double>::iterator &it);
	std::vector<double>::iterator get_rates(std::vector<double>::iterator &it);

	void print(std::ostream &out = std::cout);

};

class PSPM_Dynamic_Environment : public EnvironmentBase, public env::LightEnvironment, public env::Climate{
	public:
	
	double projected_crown_area_above_z(double t, double z, Solver *S){
		double ca_above_z = 0;
		// Loop over resident species --->
		for (int k=0; k<S->species_vec.size(); ++k){
			auto ca_above = [z,k,S](int i, double t){
				auto& p = (static_cast<Species<PSPM_Plant>*>(S->species_vec[k]))->getCohort(i);
//				double a = p.area_leaf_above(z, p.vars.height, p.vars.area_leaf);
				double ca_p = p.geometry.crown_area_extent_projected(z, p.traits);
//				std::cout << "(" << i << "," << p.geometry.get_size() << ", " << p.geometry.diameter << ", " << p.geometry.height << ", " << p.u << ", " << p.geometry.crown_area << " | " << ca_p << ")" << "\n";
				return ca_p;	
			};
//			ca_above_z += S->integrate_wudx_above(ca_above, t, z, k);
//			std::cout << "\n";
			ca_above_z += S->integrate_x(ca_above, t, k);
		}
		return ca_above_z;	
	}


	double fapar_layer(double t, int layer, Solver *S){

		double photons_abs = 0;
		for (int k=0; k<S->species_vec.size(); ++k){
			auto photons_absorbed_plant_layer = [layer,k,S, this](int i, double t){
				auto& p = (static_cast<Species<PSPM_Plant>*>(S->species_vec[k]))->getCohort(i);

				double cap_z    =            p.geometry.crown_area_above(z_star[layer], p.traits);
				double cap_ztop = (layer>0)? p.geometry.crown_area_above(z_star[layer-1], p.traits) : 0;
				double cap_layer = cap_z - cap_ztop;
				
				double Iabs_plant_layer = cap_layer * (1 - exp(-p.par.k_light * p.geometry.lai));
//				std::cout << "(" << i << "," << p.geometry.diameter << ", " << p.geometry.height << ", " << p.u << ", " << p.geometry.crown_area << ", " << cap_layer << ", " << p.geometry.lai << " | " << Iabs_plant_layer << ")" << "\n";
				return Iabs_plant_layer;
			};
//			std::cout << "\n";
			photons_abs += S->integrate_x(photons_absorbed_plant_layer, t, k);
		}		
		return photons_abs;
	}

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
	void computeEnv(double t, Solver *S){
		//updateClimate(2000+t);

		//            _xm 
		// Calculate / w(z,t)u(z,t)dz
		//        xb`
		if (use_ppa){
			z_star.clear();
			total_crown_area = projected_crown_area_above_z(t, 0, S);
			n_layers = int(total_crown_area/0.9); // Total crown projection area 

			if (n_layers < 0 || n_layers >= 50) {
				std::cout << "nlayers = " << n_layers << "\n";
				S->print(); std::cout.flush();
			}
//			if (n_layers > 5) n_layers = 5;
			assert(n_layers >= 0 && n_layers < 50);

			for (int layer = 1; layer <= n_layers; ++layer){
				auto CA_above_zstar_layer = [t, S, layer, this](double z) -> double {
					return projected_crown_area_above_z(t, z, S) - layer*0.9;
				};
				auto res = pn::zero(0, 100, CA_above_zstar_layer, 1e-4);
				z_star.push_back(res.root);
//				std::cout << "z*(" << layer << ") = " << res.root << ", CA(z*) = " << projected_crown_area_above_z(t, res.root, S) << ", " << "iter = " << res.nfnct << "\n";
			}
			z_star.push_back(0);
			//std::cout << "z*_vec (" << z_star.size() << ") = "; for(auto z: z_star) std::cout << z << " "; cout << "\n"; cout.flush();
			assert(z_star.size() == n_layers+1);
			
			canopy_openness.clear();
			fapar_tot.clear();
			canopy_openness.resize(n_layers+1);  
			fapar_tot.resize(n_layers);
			canopy_openness[0] = 1; // top layer gets 100% light
			for (int layer = 0; layer < n_layers; ++layer){  
				fapar_tot[layer] = fapar_layer(t, layer, S);
				canopy_openness[layer+1] = canopy_openness[layer] * (1-fapar_tot[layer]);
			}
			
			
		}
		else{
//			auto canopy_openness = [S, t, this](double z){
//				double kI = 0.5;

//				double leaf_area_above_z = LA_above_z(t,z,S);
//				//cout << "LA(" << z << ") = " << exp(-kI*leaf_area_above_z) << "\n";
//				return exp(-kI*leaf_area_above_z);
//			};	

//			//cout << S->xb << " " << S->getMaxSize() << endl;	
//			//time = t;
//			//for (int s=0; s<S->n_species(); ++s) S->get_species(s)->u0_save = S->get_u0(t, s);
//			light_profile.construct(canopy_openness, 0, S->maxSize());
		}	
	}
	
	
	void print(double t){
		Climate::print(t);
		LightEnvironment::print();
	}
};



#endif
