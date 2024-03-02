#include "index_utils.h"

// inline bool smaller_than(const std::vector<double>& x1, const std::vector <double>& x2){
// 	for(int i = 0; i < x1.size(); ++i){
// 		if(x1[i] >= x2[i]){
// 			return false; 
// 		}
// 	}
// 	return true;
// }

// inline bool larger_than(const std::vector<double>& x1, const std::vector <double>& x2){
// 	for(int i = 0; i < x1.size(); ++i){
// 		if(x1[i] <= x2[i]){
// 			return false; 
// 		}
// 	}
// 	return true;
// }

// inline bool larger_or_equal(const std::vector<double>& x1, const std::vector <double>& x2){
// 	for(int i = 0; i < x1.size(); ++i){
// 		if(x1[i] < x2[i]){
// 			return false; 
// 		}
// 	}
// 	return true;
// }

/// @param w            A function or function-object of the form `w(int i, double t)` that returns a double. 
///                     It can be a lambda. This function should access the 'i'th cohort and compute the weight 
///                     from the cohort's properties. The function `w` should be able access to the `Solver` in
///                     order to access cohorts. 
/// @param t            The current time (corresponding to the current physiological state). This will be passed 
///                     to `w` to allow the weights to be a direct function of time.
/// @param xlow         The lower limit of the integral
/// @param species_id   The id of the species for which the integral should be computed
///
/// \image html size_integral.png width=700cm 
///
/// the computation of this integral depends on the solver method. For different solvers, the integral is defined as follows:
/// 
/// `FMU:` \f$\quad I = \sum_{i=i_0}^J h_i w_i u_i\f$
/// 
/// `EBT:` \f$\quad I = \sum_{i=i_0}^J w_i N_i\f$, with \f$x_0 = x_b + \pi_0/N_0\f$
/// 
/// `CM :` \f$\quad I = \sum_{i=i_0}^J h_i  (w_{i+1}u_{i+1}+w_i u_i)/2\f$
/// 
/// where \f$i_0 = \text{argmax}(x_i \le x_{low})\f$, \f$h_i = x_{i+1}-x_i\f$, and \f$w_i = w(x_i)\f$. 
/// 
/// If interpolation is turned on, \f$h_{i_0}=x_{i_0+1}-x_{low}\f$, whereas \f$u(x_{low})\f$ is set 
/// to \f$u(x_{i_0})\f$ in FMU and calculated by bilinear interpolation in CM (See Figure). 
/// Interpolation does not play a role in EBT.
/// 

//             _xm 
// Calculate _/ w(z,t)u(z)dz
//         xlow
// implementation from orig plant model	
// ----
// I += (x_hi - x_lo) * (f_hi + f_lo);
// x_hi = x_lo;
// f_hi = f_lo;
// if (x_lo < xlow) break;
// ---- 
template<typename wFunc>
double Solver::integrate_wudx_above(wFunc w, double t, const std::vector<double>& xlow, int species_id){

	Species_Base* spp = species_vec[species_id];

	// cohorts are inserted at the end, hence sorted in descending order
	// FIXME: should there be an optional "sort_cohorts" control parameter? Maybe some freak models are such that cohorts dont remain in sorted order?
	// Note: CM currently only allows 1D state, so this is treated as before
	if (method == SOLVER_CM || method == SOLVER_ICM){
		// CM only handles 1D states. Cohorts must be sorted by size
		
		// integrate using trapezoidal rule 
		// Note, new cohorts are inserted at the end, so x will be descending
		bool integration_completed = false;
		double I = 0;
		double u_hi = spp->getU(0); 
		double x_hi = spp->getX(0)[0];
		double f_hi = w(0, t)*u_hi;
		
		for (int i=1; i<spp->J; ++i){  // going from high ----> low
			double u_lo = spp->getU(i); 
			double x_lo = spp->getX(i)[0];
			double f_lo = w(i, t)*u_lo;
	
			if (x_lo <= xlow[0]){  // check if xlow falls within the current interval. 
				// * if interpolation is on, stop exactly at xlow, else take the full interval
				double f = (control.integral_interpolate)? f_lo + (f_hi-f_lo)/(x_hi-x_lo)*(xlow[0] - x_lo) : f_lo;  
				double x = (control.integral_interpolate)? xlow[0] : x_lo;
				I += (x_hi - x) * (f_hi + f);
				integration_completed = true;
			}
			else{
				I += (x_hi - x_lo) * (f_hi + f_lo);
			}
		
			x_hi = x_lo;
			f_hi = f_lo;
			if (integration_completed) break;
		}
		
		//if (integration_completed) assert(x_hi < xlow);
		//if (!integration_completed) assert(x_hi >= xlow || spp.J == 1);
		
		// if integration has not completed, continue to the boundary interval	
		if (spp->J == 1 || !integration_completed){  
			// boundary at xb
			double u0 = spp->get_boundary_u();
			double x_lo = spp->xb[0];
			if (x_hi > x_lo){ 
				double f_lo =  w(-1, t)*u0; // -1 is boundary cohort
				double f = (control.integral_interpolate)? f_lo + (f_hi-f_lo)/(x_hi-x_lo)*(xlow[0] - x_lo) : f_lo;  
				double x = (control.integral_interpolate)? xlow[0] : x_lo;
				I += (x_hi-x)*(f_hi+f);
			}
		}
		
		return I*0.5;
	}

	else if (method == SOLVER_FMU || method == SOLVER_IFMU){
		//if (xlow < spp->xb) throw std::runtime_error("integral lower bound must be >= xb");
		// if (xlow > spp->x[spp->J]) return 0;  // if xlow is above the maximum size, integral is 0

		// integrate using midpoint quadrature rule
		double I=0;
		for (int i=0; i<spp->J; ++i){  
			bool x_ge_xlow = true;
			for (int k=0; k<spp->istate_size; ++k){
				x_ge_xlow = x_ge_xlow && spp->getX(i)[k] >= xlow[k];
			}
			if (x_ge_xlow){
				std::vector<double> dx = utils::tensor::coord_value(utils::tensor::index(i, spp->dim_centres), spp->h);
				double dV = std::accumulate(dx.begin(), dx.end(), 1.0, std::multiplies<double>()); // TODO: Better to precompute this and store in cohort
				I += w(i, t)*spp->getU(i)*dV;
			}
		}

		// if (xlow[0] > spp->x[0][spp->J]) return 0;  // if xlow is above the maximum size, integral is 0

		// // integrate using midpoint quadrature rule
		// for (int i=spp->J-1; i>=0; --i){  // in FMU, cohorts are sorted ascending
		// 	if (spp->x[0][i] <= xlow[0]){ // check if last interval is reached 
		// 		double h = (control.integral_interpolate)?  spp->x[0][i+1]-xlow[0]  :  spp->h[0][i];
		// 		I += h * w(i,t) * spp->getU(i);
		// 		break;
		// 	}
		// 	else{
		// 		I += spp->h[0][i] * w(i,t) * spp->getU(i);
		// 	}
		// }

		return I;
	}
	
	else if (method == SOLVER_EBT || method == SOLVER_IEBT){
		// set up cohorts to integrate
		realizeEbtBoundaryCohort(spp);

		// calculate integral
		double I = 0;
		for (int i=0; i<spp->J; ++i){  // in EBT, cohorts are sorted descending
		   	// bool x_ge_xlow = std::equal(spp->getX(i).begin(), spp->getX(i).end(),
			//                             xlow.begin(), xlow.end(),
			//                             [](int a, int b)->bool {return a >= b; });
			bool x_ge_xlow = true;
			for (int k=0; k<spp->istate_size; ++k){
				x_ge_xlow = x_ge_xlow && spp->getX(i)[k] >= xlow[k];
			}
			if (x_ge_xlow) I += w(i, t)*spp->getU(i); // if X >= xlow, we include it in the intgral
		}
		
		// restore the original pi0-cohort
		restoreEbtBoundaryCohort(spp);

		return I;
	}
	
	else if (method == SOLVER_ABM){
		// calculate integral
		double I = 0;
		for (int i=0; i<spp->J; ++i){  
			bool x_ge_xlow = true;
			for (int k=0; k<spp->istate_size; ++k){
				x_ge_xlow = x_ge_xlow && spp->getX(i)[k] >= xlow[k];
			}
			if (x_ge_xlow) I += w(i, t)*spp->getU(i); // if X >= xlow, we include it in the intgral
		}
		return I;
	}
	
	else{
		throw std::runtime_error("Unsupported solver method");
	}
}


// /// @param w            A function or function-object of the form `w(int i, double t)` that returns a double. 
// ///                     It can be a lambda. This function should access the 'i'th cohort and compute the weight 
// ///                     from the cohort's properties. The function `w` should be able access to the `Solver` in
// ///                     order to access cohorts. 
// /// @param t            The current time (corresponding to the current physiological state). This will be passed 
// ///                     to `w` to allow the weights to be a direct function of time.
// /// @param xlow         The lower limit of the integral
// /// @param species_id   The id of the species for which the integral should be computed
// ///
// /// \image html size_integral.png width=700cm 
// ///
// /// the computation of this integral depends on the solver method. For different solvers, the integral is defined as follows:
// ///
// /// `EBTN:` \f$\quad I = \sum_{i=i_0}^J w_i N_i\f$, with \f$x_0 = x_b + \pi_0/N_0\f$
// /// 
// //             _xm 
// // Calculate _/ w(z,t)u(z)dz
// //         xlow
// // implementation from orig plant model	
// // ----
// // I += (x_hi - x_lo) * (f_hi + f_lo);
// // x_hi = x_lo;
// // f_hi = f_lo;
// // if (x_lo < xlow) break;
// // ---- 
// template<typename wFunc>
// double Solver::integrate_wudxn_above(wFunc w, double t, std::vector<double> xnlow, std::vector<double> xnhigh, int species_id){

// 	Species_Base* spp = species_vec[species_id];

// 	if (method == SOLVER_EBTN || method == SOLVER_IEBTN){
// 		// set up cohorts to integrate
// 		realizeEbtnBoundaryCohort(spp);

// 		// sort cohorts, but skip pi0-cohort
// 		spp->sortCohortsDescending(1); // FIXME: Add a label to EBT boundary cohort and assert that it is always at J-1

// 		// calculate integral
// 		double I = 0;
// 		for (int i=0; i<spp->J; ++i){  // in EBT, cohorts are sorted descending
// 		   	if (smaller_than(spp->getXn(i),xnlow) && larger_than(spp->getXn(i),xnhigh)) break; // if X == xlow, we still include it in the intgral
// 			else I += w(i, t)*spp->getU(i);
// 		}
		
// 		// restore the original pi0-cohort
// 		restoreEbtnBoundaryCohort(spp);

// 		return I;
// 	}
	
// 	else{
// 		throw std::runtime_error("Unsupported solver method");
// 	}
// }





// //             _xm 
// // Calculate _/ w(z,t)u(z,t)dz
// //         xb
// // This will eventually be replaced with a call to integrate_wudx_above
// template<typename wFunc>
// double Solver::integrate_x(wFunc w, double t, int species_id){
// 	Species_Base* spp = species_vec[species_id];

// 	if (method == SOLVER_FMU || method == SOLVER_IFMU){
// 		// integrate using midpoint quadrature rule
// 		double I=0;
// 		for (unsigned int i=0; i<spp->J; ++i){
// 			I += spp->h[i]*w(i, t)*spp->getU(i);  // TODO: Replace with std::transform after profiling
// 			//std::cout << "integral: " << w(i, t) << " * " << spp->getU(i) << std::endl;
// 		}
// 		return I;
// 	}
	
// 	else if (method == SOLVER_EBT || method == SOLVER_IEBT){
// 		// integrate using EBT rule (sum over cohorts)
// 		realizeEbtBoundaryCohort(spp);

// 		// calculate integral
// 		double I = 0;
// 		for (int i=0; i<spp->J; ++i) I += w(i, t)*spp->getU(i);
		
// 		// restore the original pi0-cohort
// 		restoreEbtBoundaryCohort(spp);

// 		return I;
// 	}

// 	else if (method == SOLVER_EBTN || method == SOLVER_IEBTN){
// 		// integrate using EBT rule (sum over cohorts)
// 		realizeEbtnBoundaryCohort(spp);

// 		// calculate integral
// 		double I = 0;
// 		for (int i=0; i<spp->J; ++i) I += w(i, t)*spp->getU(i);
		
// 		// restore the original pi0-cohort
// 		restoreEbtnBoundaryCohort(spp);

// 		return I;
// 	}
	
// 	if (method == SOLVER_CM || method == SOLVER_ICM){
// 		// integrate using trapezoidal rule 
// 		// Note, new cohorts are inserted at the end, so x will be descending
// 		double I = 0;
// 		double u_hi = spp->getU(0); 
// 		double x_hi = spp->getX(0);
// 		double f_hi = w(0, t)*u_hi;

// 		for (int i=1; i<spp->J; ++i){
// 			double u_lo = spp->getU(i); 
// 			double x_lo = spp->getX(i);
// 			double f_lo = w(i, t)*u_lo;
	
// 			I += (x_hi - x_lo) * (f_hi + f_lo);
// 			x_hi = x_lo;
// 			f_hi = f_lo;
// 		}
		
// 		// boundary at xb
// 		double u0 = spp->get_boundary_u();
// 		double x_lo = spp->xb;
// 		double f_lo =  w(-1, t)*u0; // -1 is boundary cohort
// 		I += (x_hi-x_lo)*(f_hi+f_lo);
		
// 		return I*0.5;
// 	}
	
// 	else if (method == SOLVER_ABM){
// 		// calculate integral. Sorting is not required here because all cohorts will be touched anyway
// 		double I = 0;
// 		for (int i=0; i<spp->J; ++i){       // in ABM, cohorts are sorted descending
// 			I += w(i, t)*spp->getU(i);
// 		}
		
// 		return I;
// 	}
	
// 	else{
// 		throw std::runtime_error("Unsupported solver method");
// 	}
// }


//             _xm 
// Calculate _/ w(z,t)u(z,t)dz
//         xb
// This will eventually be replaced with a call to integrate_wudx_above
template<typename wFunc>
double Solver::state_integral(wFunc w, double t, int species_id){
	Species_Base* spp = species_vec[species_id];

	// NOTE: This works for 1D state only! Cohorts assumed to be sorted
		if (method == SOLVER_CM || method == SOLVER_ICM){
		// spp->sortCohortsDescending(0);
		// copyCohortsToState();
		// integrate using trapezoidal rule 
		// Note, new cohorts are inserted at the end, so x will be descending
		double I = 0;
		double u_hi = spp->getU(0); 
		double x_hi = spp->getX(0)[0];
		double f_hi = w(0, t)*u_hi;
		// std::cout << "size_integral: x = " << x_hi << ", w(0,t) = " << w(0,t) << '\n';
		for (int i=1; i<spp->J; ++i){
			double u_lo = spp->getU(i); 
			double x_lo = spp->getX(i)[0];
			double f_lo = w(i, t)*u_lo;
			if (debug) std::cout << "size_integral: x = " << x_lo << ", w(" << i << ",t) = " << w(i,t) << '\n';
	
			I += (x_hi - x_lo) * (f_hi + f_lo);
			x_hi = x_lo;
			f_hi = f_lo;
		}
		
		// boundary at xb
		double u0 = spp->get_boundary_u();
		double x_lo = spp->xb[0];
		double f_lo =  w(-1, t)*u0; // -1 is boundary cohort
		I += (x_hi-x_lo)*(f_hi+f_lo);
		
		return I*0.5;
	}

	if (method == SOLVER_EBT || method == SOLVER_IEBT){
		// integrate using EBT rule (sum over cohorts)
		realizeEbtBoundaryCohort(spp);

		// calculate integral
		double I = 0;
		for (int i=0; i<spp->J; ++i){
			// auto x = spp->getX(i);
			// if (larger_or_equal(x, xnlow) && smaller_than(x, xnhigh)){
				I += w(i, t)*spp->getU(i);
			// }
		}
		
		// restore the original pi0-cohort
		restoreEbtBoundaryCohort(spp);

		return I;
	}

	if (method == SOLVER_FMU || method == SOLVER_IFMU){
		// integrate using midpoint quadrature rule
		double I=0;
		for (unsigned int i=0; i<spp->J; ++i){
			std::vector<double> dx = utils::tensor::coord_value(utils::tensor::index(i, spp->dim_centres), spp->h);
			double dV = std::accumulate(dx.begin(), dx.end(), 1.0, std::multiplies<double>());

			I += w(i, t)*spp->getU(i)*dV;  // TODO: Replace with std::transform after profiling
			//std::cout << "integral: " << w(i, t) << " * " << spp->getU(i) << std::endl;
		}
		return I;
	}

	else if (method == SOLVER_ABM){
		// calculate integral
		double I = 0;
		for (int i=0; i<spp->J; ++i){  // in EBT, cohorts are sorted descending
			I += w(i, t)*spp->getU(i); // if X >= xlow, we include it in the intgral
		}
		return I;
	}

	else{
		throw std::runtime_error("Unsupported solver method");
	}
}
