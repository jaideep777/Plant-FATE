#include <iomanip>
#include "pspm_dynamic_environment.h"
#include "pspm_interface.h"
#include "trait_evolution.h"

using namespace std;



// **********************************************************
// ********** PSPM_Dynamic_Environment **********************
// **********************************************************


PSPM_Dynamic_Environment::PSPM_Dynamic_Environment() : LightEnvironment() {
	
}

/// @ingroup    ppa_module
/// @brief      Calculate the total crown area above height z, contributed by all individuals of the resident species.
/// @param t    Time in current timestep.
/// @param z    distance from the ground
/// @param S    Pointer to the Solver being used.
/// @return     Total crown area above z
/// @details    The crown area above z from all individuals of species `k` is calculated as \f[A_k = \int_{x_b}^{x_m}{A_{cp}(s)u(s)ds},\f]
///             where \f$A_{cp}\f$ is the projected crown area at height z, including area covered by gaps. This is
///             calculated by the function plant::PlantGeometry::crown_area_extent_projected.  
///             The total crown area above z is the \f[A = \sum_k {A_k}\f]
double PSPM_Dynamic_Environment::projected_crown_area_above_z(double t, double z, Solver *S){
	double ca_above_z = 0;
	// Loop over resident species --->
	for (int k=0; k<S->species_vec.size(); ++k){

		// skip mutants
		auto spp = static_cast<MySpecies<PSPM_Plant>*>(S->species_vec[k]);
		if (!spp->isResident) continue;

		auto ca_above = [z,spp](int i, double t){
			auto& p = spp->getCohort(i);
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

/// @ingroup    ppa_module
double PSPM_Dynamic_Environment::fapar_layer(double t, int layer, Solver *S){

	double photons_abs = 0;
	for (int k=0; k<S->species_vec.size(); ++k){
		
		// skip mutants
		auto spp = static_cast<MySpecies<PSPM_Plant>*>(S->species_vec[k]);
		if (!spp->isResident) continue;

		auto photons_absorbed_plant_layer = [layer,spp, this](int i, double t){
			auto& p = spp->getCohort(i);

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


/// @brief        Solver interface for updating the environment from the given state.
/// @param t      Time in current timestep.
/// @param S      Pointer to the Solver being used, provided for computing state integrals. The Solver provides this via a `this` reference. 
/// @param _S     Iterator to the `state` vector used by the ODE solver.
/// @param _dSdt  Iterator to the `rates` vector used by the ODE solver.
/// @ingroup      ppa_module
/// @details 
/// This function must perform a complete update of the Environment, 
/// including any species or individual specific differences in the env. 
/// If Env is a function of E(x), the solver's SubdivisionSpline can be used to store E(x).
///
/// This function receives the updated state vector from the Solver's ODE stepper,
/// from which the complete system state can be retrieved.
/// This is the only function to which the Solver exposes the state vector, so if desired, the state vector 
/// can be saved from here and reused in other rate functions.
/// 
/// Interface to the rates vector can be used for computing the rates of ODE-based environmental 
/// variables if any such have been added as system variables - e.g., species-level seed pools can be implemented 
/// through this mechanism. 
void PSPM_Dynamic_Environment::computeEnv(double t, Solver *S, std::vector<double>::iterator _S, std::vector<double>::iterator _dSdt){
	updateClimate(t);

	//            _xm 
	// Calculate / w(z,t)u(z,t)dz
	//        xb`
	if (use_ppa){
		double fG = 0.99;
		z_star.clear();
		total_crown_area = projected_crown_area_above_z(t, 0, S);
		n_layers = int(total_crown_area/fG); // Total crown projection area 

		if (n_layers < 0 || n_layers >= 50) {
			std::cout << "nlayers = " << n_layers << "\n";
			S->print(); std::cout.flush();
		}
//			if (n_layers > 5) n_layers = 5;
		assert(n_layers >= 0 && n_layers < 50);

		for (int layer = 1; layer <= n_layers; ++layer){
			auto CA_above_zstar_layer = [t, S, layer, fG, this](double z) -> double {
				return projected_crown_area_above_z(t, z, S) - layer*fG;
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
		throw std::runtime_error("Only PPA mode is implemented currently. Set use_ppa to true");
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


void PSPM_Dynamic_Environment::print(double t){
	Climate::print(t);
	LightEnvironment::print();
}
