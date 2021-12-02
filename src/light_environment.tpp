namespace env{

template <class Community>
double LightEnvironment::projected_crown_area_above_z(double t, double z, Community *C){
	double ca_above_z = 0;
	//// Loop over resident species --->
	//for (int k=0; k<S->species_vec.size(); ++k){
		//auto la_above = [z,k,S](int i, double t){
			//auto& p = ((Species<PSPM_Plant>*)S->species_vec[k])->getCohort(i);
			//double a = p.area_leaf_above(z, p.vars.height, p.vars.area_leaf);
			////cout << "(" << i << "," << p.u << ", " << a << ")" << "\t";
			//return a;	
		//};
		//leaf_area_above_z += S->integrate_wudx_above(la_above, t, z, k);
		////cout << "\n";
		////leaf_area_above_z += S->integrate_x(la_above, t, state_vec, i);
	//}
	return ca_above_z;
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
template <class Community>
void LightEnvironment::computeEnv(double t, Community *C){
//	//            _xm 
//	// Calculate / w(z,t)u(z,t)dz
//	//        xb`
//	if (use_ppa){
//		z_star.clear();
//		n_layers = int(LA_above_z(t, 0, S)); // Get total leaf area, = area above height 0.
//		if (n_layers < 0 || n_layers >= 50) {
//			std::cout << "nlayers = " << n_layers << "\n";
//			S->print();
//		}
//		assert(n_layers >= 0 && n_layers < 50);
//		for (int layer = 1; layer <= n_layers; ++layer){
//			auto LA_l = [t, S, layer, this](double z) -> double {
//				return LA_above_z(t, z, S)-layer; 
//			};
//			auto res = pn::zero(0, S->maxSize(), LA_l, 1e-4);
//			z_star.push_back(res.root);
//			//cout << "z*(" << layer << ") = " << res.root << ", " << "iter = " << res.nfnct << "\n";
//		}
//		z_star.push_back(0);
//		//std::cout << "z*_vec (" << z_star.size() << ") = "; for(auto z: z_star) std::cout << z << " "; cout << "\n"; cout.flush();
//		assert(z_star.size() == n_layers+1);
//	}
//	else{
//	auto canopy_openness = [S, t, this](double z){
//		double kI = 0.5;

//		double leaf_area_above_z = LA_above_z(t,z,S);
//		//cout << "LA(" << z << ") = " << exp(-kI*leaf_area_above_z) << "\n";
//		return exp(-kI*leaf_area_above_z);
//	};	

//	//cout << S->xb << " " << S->getMaxSize() << endl;	
//	//time = t;
//	//for (int s=0; s<S->n_species(); ++s) S->get_species(s)->u0_save = S->get_u0(t, s);
//	light_profile.construct(canopy_openness, 0, S->maxSize());
//	}

}

} // namespace env


