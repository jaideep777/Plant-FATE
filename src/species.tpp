#include <cassert>
#include <algorithm>
#include <numeric>
#include <iomanip>

#include "io_utils.h"
#include "fof.h"

// *************** Species_Base   ***************  

template<class T>
void Species_Base::addCohort(T bc){
	(dynamic_cast<Species<T>>(*this)).addCohort(bc);
}


// *************** Species<Model> ******************

// Actually, this constructor is not required because it is covered by the second one below
// copy-construct boundary cohort
// template<class Model>
// Species<Model>::Species(const Model& M) : boundaryCohort(M){
// 	istate_size = M.state_size;
// }

// construct boundary cohort using default / custom / copy constructor
template<class Model>
template <typename... ARGS>
Species<Model>::Species(ARGS... args) : boundaryCohort(args...){
	istate_size = boundaryCohort.state_size;
}


template<class Model>
void Species<Model>::clear_vectors(){
	dim_centres.clear();
	dim_edges.clear();
	X.clear();
	x.clear();
	h.clear();
	xb.clear();
	cohorts.clear();
}


template<class Model>
void Species<Model>::resize(int _J){
	J = _J;
	cohorts.resize(J, boundaryCohort);  // when resizing always insert copies of boundary cohort
}

// TODO: make this comment more sensible
// Returns an "artificial" vector of largest xk in all the cohorts - no cohort has to actually occupy this point though 
// We dont know in advance whether there is a boundary cohort at position J-1. So provide an option to skip it
template<class Model>
std::vector<double> Species<Model>::get_maxSize(int skip){ // TODO ALERT: make sure this sees the latest state
	if (cohorts.empty()) return xb;

	std::vector<double> largest = xb;
	for (size_t i = 0; i < J-skip; ++i){ 
		for(size_t k = 0; k < largest.size(); ++k){
			largest[k] = std::max(cohorts[i].x[k], largest[k]);
		}
	}
	return largest;	
}


// // FIXME: this constructor should  be corrected or removed.
// template<class Model>
// Species<Model>::Species(std::vector<double> breaks){
// 	istate_size = Model::state_size;
// 	J = breaks.size(); // changes with 
// 	cohorts.resize(J, boundaryCohort);
// 	for (int i=0; i<J; ++i) cohorts[i].x = breaks[i];	
// }




template <class Model>
void Species<Model>::print(){
	std::cout << "~~~~~ Species ~~~~~\n";
	std::cout << "xb = " << boundaryCohort.x << " / " << xb << "\n";
	std::cout << "xsize = " << J << "\n";
	std::cout << "istate size = " << istate_size << '\n';
	std::cout << "Extra state variables: " << n_accumulators << "\n";
	std::cout << "Total state size: " << stateSizeTotal() << "\n";
	std::cout << "Input birth flux = " << birth_flux_in << "\n";

	print_extra();

	std::cout << "x (" << x.size() << "):\n"; 
	for (auto xx : x) std::cout << xx << "\n";
	std::cout << "X (" << X.size() << "):\n"; 
	for (auto xx : X) std::cout << xx << "\n";
	std::cout << "h (" << h.size() << "):\n"; 
	for (auto xx : h) std::cout << xx << "\n";

	std::cout << "Cohorts: (" << cohorts.size() << ")\n";
	std::cout << std::setw(6) << "t0";
	for (int i=0; i<istate_size; ++i){
		std::cout << std::setw(12) << "x[" << i << "]"; 
	}
	std::cout << std::setw(12) << "u" << " ";
	if (!cohorts.empty()){
		for (auto s : cohorts[0].varnames) std::cout << std::setw(12) << s << " ";
	}
	std::cout << "\n";
	for (auto& c : cohorts) {
		//c.print_xu(); //std::cout << c.x << "\t" << c.u << "\t";
		c.print();
		std::cout << "\n";
	}
	std::cout << "- - - - - - - - - - - - - - - - - - - - - - - \n";
	//boundaryCohort.print_xu(); //std::cout << c.x << "\t" << c.u << "\t";
	boundaryCohort.print();
	std::cout << "\n";
	std::cout << "Max size = " << get_maxSize() << "\n";
	std::cout << "-------\n\n"; std::cout.flush();

}


template <class Model>
Cohort<Model>& Species<Model>::getCohort(int i){
	if (i == -1) return boundaryCohort;
	else return cohorts[i];
}

template <class Model>
std::vector<double> Species<Model>::getX(int i){
	if (i == -1) return to_vector(boundaryCohort.x);
	else return to_vector(cohorts[i].x);
}


template <class Model>
double Species<Model>::getU(int i){
	return cohorts[i].u;
}


template <class Model>
double Species<Model>::get_boundary_u(){
	return boundaryCohort.u;
}


template <class Model>
void Species<Model>::set_xb(std::vector<double> _xb){
	xb = _xb;
	boundaryCohort.set_size(_xb);
}

template <class Model>
void Species<Model>::set_ub(double _ub){
	boundaryCohort.u = _ub;
}


template <class Model>
void Species<Model>::set_birthTime(int i, double t0){
	cohorts[i].birth_time = t0;
}


template <class Model>
void Species<Model>::setX(int i, std::vector<double> _x){
	cohorts[i].set_size(_x);
}

template <class Model>
void Species<Model>::setU(int i, double _u){
	cohorts[i].u = _u;
}

// template <class Model>
// std::vector <double> Species<Model>::getStateAt(int i){
// 	std::vector<double> _xn = cohorts[i].xn;
// 	return(_xn);
// }

// //To remove
// template <class Model>
// double Species<Model>::dXn(int i){
// 	double dxn = 1;
// 	for (size_t k = 0; k < xn[i].size(); ++k){
// 		double dx_k = (xn[i][k] - next_xk_desc(xn[i][k], k));
// 		dxn = dxn * dx_k;
// 	}
// 	return dxn;
// }

// template <class Model> 
// double Species<Model>::dXn(std::vector<double> xn1, std::vector<double> xn2){
// 	double _dxn = 1;
// 	for (size_t k = 0; k < xn1.size(); ++k){
// 		double dx_k = (xn2[k] - xn1[k]);
// 		_dxn = _dxn * dx_k;
// 	}
// 	return _dxn;
// }

// template <class Model> 
// std::vector<double> Species<Model>::cohort_dist(std::vector<double> xn1, std::vector<double> xn2){
// 	std::vector<double> _dxn;
// 	for (size_t k = 0; k < xn1.size(); ++k){
// 		double dx_k = (xn2[k] - xn1[k]);
// 		_dxn.push_back(dx_k);
// 	}
// 	return _dxn;
// }

// template <class Model>
// double Species<Model>::next_xk_desc(double xnk, int k){
// 	double next_smallest = xnb[k];
// 	for(size_t i = 0; i < xn.size(); i++){
// 		if(xn[i][k] < xnk && xn[i][k] > next_smallest){
// 			next_smallest = xn[i][k];
// 		}
// 	}
// 	return next_smallest;
// }

// template <class Model>
// double Species<Model>::next_xk_asc(double xnk, int k){
// 	double next_biggest = get_maxSizeN()[k];
// 	for(size_t i = 0; i < xn.size(); i++){
// 		if(xn[i][k] > xnk && xn[i][k] < next_biggest){
// 			next_biggest = xn[i][k];
// 		}
// 	}
// 	return next_biggest;
// }

// // FIXME: need to rethink this... 
// template <class Model>
// std::vector<double> Species<Model>::next_xn_desc(std::vector<double> _xni){
// 	std::vector<double> next_smallest = xnb; //boundary is always smaller
// 	for (size_t i = 1; i <= J; ++i){ //should I check the boundary cohort too? probably not needed can be i < J
// 		bool smaller = true;
// 		for(size_t k = 0; k < next_smallest.size(); ++k){
// 				if(xn[i][k] < _xni[k] && xn[i][k] > next_smallest[k]){
// 					smaller = false; 
// 					break;
// 				}
// 		}
// 		if(smaller){
// 			next_smallest = cohorts[i].xn;
// 		}
// 	}
// 	return next_smallest;
// }

// // FIXME: need to rethink this... 
// template <class Model>
// std::vector<double> Species<Model>::next_xn_asc(std::vector<double> _xni){
// 	std::vector<double> next_biggest = get_maxSizeN();
// 	for (size_t i = 1; i <= J; ++i){ //should I check the boundary cohort too? probably not needed can be i < J
// 		bool bigger = true;
// 		for(size_t k = 0; k < next_biggest.size(); ++k){
// 			if(xn[i][k] > _xni[k] && xn[i][k] < next_biggest[k]){
// 				bigger = false; 
// 				break;
// 			}
// 		}
// 		if(bigger){
// 			next_biggest = cohorts[i].xn;
// 		}
// 	}
// 	return next_biggest;
// }


// ---------------------------

template <class Model>
void Species<Model>::initAccumulators(double t, void * env){
	// init boundary cohort 
	boundaryCohort.init_accumulators(t, env); 
	// init internal cohorts
	for (auto& c : cohorts){
		c.init_accumulators(t, env);	// init state
	}
}


template <class Model>
void Species<Model>::initAndCopyAccumulators(double t, void * env, std::vector<double>::iterator &it){
	// init boundary cohort (no copy required)
	boundaryCohort.init_accumulators(t, env); 
	// init internal cohorts and copy to state vector
	for (auto& c : cohorts){
		c.init_accumulators(t, env);	// init state
		auto it_prev = it;		
		it = c.get_accumulators(it);	// copy the initialized state into state vector
		assert(distance(it_prev, it) == n_accumulators);
	}
}


template <class Model>
void Species<Model>::initBoundaryCohort(double t, void * env){
	boundaryCohort.birth_time = t;
	boundaryCohort.init_accumulators(t, env);
}


template <class Model>
double Species<Model>::init_density(int i, void * env){
	assert(i>=0);
	return cohorts[i].init_density(env, birth_flux_in);
}

// TODO: check increment here itself
template <class Model>
void Species<Model>::copyAccumulatorsToCohorts(std::vector<double>::iterator &it){
	for (auto& c : cohorts) c.set_accumulators(it);
}


template <class Model>
void Species<Model>::copyAccumulatorsToState(std::vector<double>::iterator &it){
	for (auto& c : cohorts) c.get_accumulators(it);
}

template <class Model>
void Species<Model>::accumulatorRates(std::vector<double>::iterator &it){
	for (auto& c : cohorts) c.get_accumulatorRates(it);
}


template <class Model>
void Species<Model>::triggerPreCompute(){
	for (auto& c : cohorts) c.need_precompute = true;
	boundaryCohort.need_precompute = true;
}


template <class Model>
double Species<Model>::establishmentProbability(double t, void * env){
	return boundaryCohort.establishmentProbability(t, env);
}


template <class Model>
double Species<Model>::calc_boundary_u(std::vector<double> gb, double pe){
	//std::cout << "calc_boundary_u\n";
	if (bfin_is_u0in){
		boundaryCohort.u = birth_flux_in;
	}
	else {
		double gdotn = std::accumulate(gb.begin(), gb.end(), 0);  // FIXME JJ: Check if this is correct
		boundaryCohort.u = (gdotn>0)? birth_flux_in * pe/gdotn  :  0;  
	}
	return boundaryCohort.u;
}


template <class Model>
std::vector<double> Species<Model>::growthRate(int i, double t, void * env){
	Cohort<Model> &c = (i<0)? boundaryCohort : cohorts[i];
	return to_vector(c.growthRate(t,env));
}


template <class Model>
std::vector<double> Species<Model>::growthRateOffset(int i, const std::vector<double>& x, double t, void * env){
	Cohort<Model> coff = (i<0)? boundaryCohort : cohorts[i];
	coff.set_size(x);
	//coff.preCompute(coff.x,t,env);

	return to_vector(coff.growthRate(t,env));
}


// Return format:
//    [ k--->
//  j   [ g1       g2       gk       g4       ... ]
//  |   [ --------------------------------------- ]
//  |   [ dg1_dx1  dg2_dx1  ...      dg4_dx1  ... ]
//  v   [ dg1_dx2  dg2_dx2  dgk_dxj  dg4_dx2  ... ]
//      [ ...                                     ]
//    ]
template <class Model>
std::vector<std::vector<double>> Species<Model>::growthRateGradient(int i, double t, void * env, const std::vector<double>& grad_dx){

	Cohort<Model> &c = (i<0)? boundaryCohort : cohorts[i];

	std::vector<std::vector<double>> g_gx;
	g_gx.reserve(istate_size+1);

	std::vector<double> g = to_vector(c.growthRate(t,env));
	g_gx.push_back(g);

	for(int k=0; k < istate_size; ++k){
		Cohort<Model> cplus = c;
		std::vector<double> xplus(c.x.begin(), c.x.end());
		xplus[k] += grad_dx[k];	// FIXME: JJ: should we normalize here by maxSize[k]?
		cplus.set_size(xplus);
		std::vector<double> gplus_k = to_vector(cplus.growthRate(t, env));
		std::vector<double> gx(istate_size);
		for (int j=0; j<istate_size; ++j){
			gx[j] = (gplus_k[j] - g[j])/grad_dx[k];
		}

		g_gx.push_back(gx);
	}
	
	return g_gx;
}


// // template <class Model>
// // std::vector<double> Species<Model>::growthRateGradientCentered(int i, double xplus, double xminus, double t, void * env){
// // 	assert(i>=0);

// // 	Cohort<Model> cplus = cohorts[i];
// // 	cplus.set_size(xplus);
// // 	//cplus.preCompute(cplus.x,t,env);

// // 	Cohort<Model> cminus = cohorts[i];
// // 	cminus.set_size(xminus);
// // 	//cminus.preCompute(cminus.x,t,env);
	
// // 	double gplus  = cplus.growthRate(cplus.x, t, env);
// // 	double gminus = cminus.growthRate(cminus.x, t, env);
	
// // 	return {gplus, gminus};
// // }


template <class Model>
double Species<Model>::mortalityRate(int i, double t, void * env){
	assert(i>=0);
	Cohort<Model> &c = cohorts[i];
	return c.mortalityRate(t,env);
}

// Return format:
//    [
//     m
//     ---
//     dm_dx1
//     dm_dx2
//     ...
//    ]
template <class Model>
std::vector<double> Species<Model>::mortalityRateGradient(int i, double t, void * env, const std::vector<double>& grad_dx){
	Cohort<Model> &c = (i<0)? boundaryCohort : cohorts[i];

	std::vector<double> m_dmdx;
	m_dmdx.reserve(istate_size+1);

	double m = c.mortalityRate(t,env);
	m_dmdx.push_back(m);

	for(int k=0; k < istate_size; ++k){
		Cohort<Model> cplus = c;
		std::vector<double> xplus(c.x.begin(), c.x.end());
		xplus[k] += grad_dx[k];	// FIXME: JJ: should we normalize here by maxSize[k]?
		cplus.set_size(xplus);
		double mplus_k = cplus.mortalityRate(t, env);
		m_dmdx.push_back( (mplus_k-m)/grad_dx[k] );
	}
	
	return m_dmdx;
}
	

template <class Model>
double Species<Model>::birthRate(int i, double t, void * env){
	Cohort<Model> &c = (i<0)? boundaryCohort : cohorts[i];
	return c.birthRate(t,env);
}	



template <class Model>
void Species<Model>::addCohort(int n){
	cohorts.reserve(cohorts.size()+n);
	for (int i=0; i<n; ++i){
		cohorts.push_back(boundaryCohort);
		++J;
	}
}


// This function allows a user to add a custom cohort to the species any time.
// However, Solver->resizeStateFromSpecies() must be called after calling this function.
template <class Model>
template <class T>
void Species<Model>::addCohort(T bc){
	cohorts.push_back(bc);
	++J;
	// FIXME: add option to sort cohorts here.
}


template <class Model>
void Species<Model>::markCohortForRemoval(int i){
	cohorts[i].remove = true;
}

template <class Model>
void Species<Model>::removeMarkedCohorts(){
	// remove marked cohorts
	auto it_end = std::remove_if(cohorts.begin(), cohorts.end(), [](Cohort<Model> &c){return c.remove;});
	cohorts.erase(it_end, cohorts.end());
	
	// reset size
	J = cohorts.size();
}


template <class Model>
void Species<Model>::markDensestCohort(){
//	group_cohorts(cohorts, 1e-6);
// 	if (cohorts.size() < 3) return; // do nothing if there are 2 or fewer cohorts
// 	int i_min = 0;
// 	double dx_min = dXn(next_xn_asc(cohorts[i_min].xn), next_xn_desc(cohorts[i_min].xn));
// 	std::vector<double> maxCohort = get_maxSizeN();
// 	if(cohorts[i_min].xn == maxCohort){ //this should be executed at most twice
// 		i_min = ++i_min;
// 		double dx_min = dXn(next_xn_asc(cohorts[i_min].xn), next_xn_desc(cohorts[i_min].xn));
// 	}

// 	for (int i=(i_min + 1); i<J-1; ++i){ // skip first and last cohorts
// 		if(maxCohort == cohorts[i].xn){
// 			continue;
// 		}
// 		double dx = dXn(next_xn_asc(cohorts[i].xn), next_xn_desc(cohorts[i].xn));
// 		if (dx < dx_min){
// 			dx_min = dx;
// 			i_min = i;
// 		}
// 	}

// 	//std::cout << "Removing cohort no. " << i_min << std::endl;
// 	cohorts.erase(cohorts.begin()+i_min);
// 	--J;
}


// Jaideep FIXME: Maybe dxcut need not be a vector... in any case, it will be hard to specify it as a control parameter when dim is not known. 
// Instead, we can use a relative dxcut across all dims
template <class Model>
void Species<Model>::markDenseCohorts(double dxcut){
	// group cohorts
	group_cohorts(cohorts, dxcut);
	
	// sort by group_id, and within group_id, sort ascending by birth_time
	std::sort(cohorts.begin(), cohorts.end(), 
		[](const Cohort<Model>&c1, const Cohort<Model>& c2){ 
			return (c1.group_id  < c2.group_id) || 
			       (c1.group_id == c2.group_id && c1.birth_time < c2.birth_time);
		}
	);

	// Mark 1st (oldest) cohort within each group to remove, skip 1-sized groups
	int prev_group_id = -999;
	for (int i = 0; i < cohorts.size(); ++i) {
		if (cohorts[i].group_id != prev_group_id) {
			if (cohorts[i].group_size > 1) {
				cohorts[i].remove = true; // Mark the first Cohort in this group
			}
			prev_group_id = cohorts[i].group_id; // Update the previous group_id
		}
    }	

	// // mark cohorts to remove; skip 1st and last cohort
	// if (cohorts.size() < 3) return; // do nothing if there are 2 or fewer cohorts
	// std::vector<double> maxCohort = get_maxSizeN();
	
	// for (int i=0; i<J-1; i+=2){
	// 	if(maxCohort == cohorts[i].x){
	// 		continue;
	// 	}
	// 	std::vector<double> dx_lo = cohort_dist(next_xn_asc(cohorts[i].xn), cohorts[i].xn);
	// 	std::vector<double> dx_hi = cohort_dist(cohorts[i].xn, next_xn_desc(cohorts[i].xn));

	// 	if (dx_lo < dxcut || dx_hi < dxcut) cohorts[i].remove = true;
	// }

	// // remove marked cohorts
	// auto it_end = std::remove_if(cohorts.begin(), cohorts.end(), [](Cohort<Model> &c){return c.remove;});
	// cohorts.erase(it_end, cohorts.end());

	// // reset size
	// J = cohorts.size();
}

template <class Model>
void Species<Model>::markDeadCohorts(double ucut){
	// mark cohorts to remove; skip pi0-cohort (index J-1)
	for (int i=0; i<J-1; ++i){
		if (cohorts[i].u < ucut) cohorts[i].remove = true;
	}
}


// // Jaideep FIXME: Maybe dxcut need not be a vector... in any case, it will be hard to specify it as a control parameter when dim is not known. 
// // Instead, we can use a relative dxcut across all dims
// template <class Model>
// void Species<Model>::mergeCohortsAddU(std::vector<double> dxcut){
// 	// // mark cohorts to remove; skip 1st and last cohort
// 	// std::vector<double> maxCohort = get_maxSizeN();
// 	// for (int i=0; i<J-1; i+=2){
// 	// 	if(maxCohort == cohorts[i].xn){
// 	// 		continue;
// 	// 	}
// 	// 	std::vector<double> dx = cohort_dist(cohorts[i].xn, next_xn_asc(cohorts[i].xn));
// 	// 	if (dx < dxcut){
// 	// 		cohorts[i-1].remove = true;
// 	// 		// FIXME: Need to also average extra state?
// 	// 		for(int k=0; k<xnb.size(); ++k){
// 	// 			cohorts[i].xn[k] = (cohorts[i].xn[k]*cohorts[i].u + cohorts[i-1].xn[k]*cohorts[i-1].u)/(cohorts[i].u + cohorts[i-1].u);
// 	// 		}

// 	// 		cohorts[i].u = cohorts[i].u + cohorts[i-1].u;
// 	// 	}
// 	// }

// 	// // remove marked cohorts
// 	// auto it_end = std::remove_if(cohorts.begin(), cohorts.end(), [](Cohort<Model> &c){return c.remove;});
// 	// cohorts.erase(it_end, cohorts.end());

// 	// // reset size
// 	// J = cohorts.size();

// }


// in EBT, one can set skip=1 to exclude the boundary cohort (J-1) from sorting
// It's best to always keep it in position J-1, otherwise, if it is realized and 
// then cohorts are sorted, it will be hard to keep track of where it went.
template <class Model>
void Species<Model>::sortCohortsDescending(size_t dim, int skip){
	assert(dim <= istate_size);
	std::sort(cohorts.begin(), cohorts.end()-skip, [dim](const Cohort<Model> &a, const Cohort<Model> &b){return a.x[dim] > b.x[dim];});
}

template <class Model>
void Species<Model>::sortCohortsAscending(size_t dim, int skip){
	assert(dim <= istate_size);
	std::sort(cohorts.begin(), cohorts.end()-skip, [dim](const Cohort<Model> &a, const Cohort<Model> &b){return a.x[dim] < b.x[dim];});
}


template <class Model>
void Species<Model>::save(std::ostream &fout){
//	Species_Base::save(fout);
	fout << "Species<T>::v2\n";
	fout << std::make_tuple(
		J
	  , istate_size
	  , n_accumulators
	  , noff_abm
	  , birth_flux_in
	  , bfin_is_u0in);
	fout << '\n';
	fout << xb << '\n';
	fout << X.size() << ' ' << x.size() << ' ' << h.size() << '\n';
	for (int i=0; i<X.size(); ++i) fout << X[i] << '\n';
	for (int i=0; i<x.size(); ++i) fout << x[i] << '\n';
	for (int i=0; i<h.size(); ++i) fout << h[i] << '\n';

	boundaryCohort.save(fout, n_accumulators);
	for (auto& C : cohorts) C.save(fout, n_accumulators);
}

template <class Model>
void Species<Model>::restore(std::istream &fin){
//	Species_Base::restore(fin);
	std::cout << "Restoring Species<T>" << std::endl;
	std::string s; fin >> s; // version number (discard)
	assert(s == "Species<T>::v2");
	fin >> J
	    >> istate_size
	    >> n_accumulators
	    >> noff_abm
	    >> birth_flux_in
	    >> bfin_is_u0in;
	
	fin >> xb;

	int nx;
	fin >> nx; X.resize(nx);
	fin >> nx; x.resize(nx);
	fin >> nx; h.resize(nx);
	for (int i=0; i<X.size(); ++i) fin >> X[i];
	for (int i=0; i<x.size(); ++i) fin >> x[i];
	for (int i=0; i<h.size(); ++i) fin >> h[i];

	boundaryCohort.restore(fin, n_accumulators);
	cohorts.resize(J, boundaryCohort); // cohorts must always be copy-constructed from the boundary cohort

	for (auto& C : cohorts) C.restore(fin, n_accumulators);
}


//TODO: maybe fix this later
template <class Model>
void Species<Model>::printCohortVector(std::ostream &out){

	for (int i=0; i<istate_size; ++i) out << "x[" << i << "]" << '\t';
	out << "u\n";
	// std::cout << "J is " << J << std::endl;
	for(int i=0; i < cohorts.size(); ++i){
		for (auto xx : cohorts[i].x) out << xx << "\t";
		out << cohorts[i].u << "\n";
	}
}

// //template <class Model>
// //void Species<Model>::backupCohort(int j){
// //	savedCohort = cohorts[j];
// //}

// //template <class Model>
// //void Species<Model>::restoreCohort(int j){
// //	cohorts[j] = savedCohort;
// //}

// //template <class Model>
// //void Species<Model>::copyBoundaryCohortTo(int j){
// //	cohorts[j] = boundaryCohort;
// //}

// //template <class Model>
// //void Species<Model>::backupBoundaryCohort(){
// 	//boundaryCohort_backup = boundaryCohort;
// //}

// //template <class Model>
// //void Species<Model>::restoreBoundaryCohort(){
// 	//boundaryCohort = boundaryCohort_backup;
// //}


