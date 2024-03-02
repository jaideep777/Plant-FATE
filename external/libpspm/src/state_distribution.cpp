#include "solver.h"
#include "index_utils.h"
#include <algorithm>


// xb
// |---*--|--*--|----*----|---*---|
//     0     1       2        3        <--- points
// 0      1     2         3       4    <--- breaks

// Test this function with this R script:
/**
x = c(0, 1,1.1,1.2, 2,2.2,2.3, 4.5,5.5,6.5)[length(x):1]
N = c(1, 1,1,  1,   2,2,  2,   3,   4,  5)[length(N):1]
plot(N~x)
abline(v=breaks, col="grey")
breaks = seq(0,10,1)
dens = rep(0, length(x))
J = length(x)

current_interval = length(breaks)-1
for (i in 1:J){
  while(breaks[current_interval] > x[i]) current_interval = current_interval-1
  dens[current_interval] = dens[current_interval]+N[i]
}
**/
struct point{
	double xmean = 0;
	double abund = 0;
	int    count = 0;
	double h = 0;
	std::vector<double> xmean_vec;
};	

// breaks must be sorted ascending
std::vector<double> Solver::getDensitySpecies1D(int k, int dim, const std::vector<double>& breaks, Spline::Extr extrapolation_method){
	auto spp = species_vec[k];
	std::vector <double> xx, uu;
	
	// if (method == SOLVER_ABM) spp->sortCohortsDescending(dim); // ABM needs sorted cohorts
	// if (method == SOLVER_EBT || method == SOLVER_IEBT) spp->sortCohortsDescending(dim, 1); // for EBT, skip boundary cohort and sort

	if (method == SOLVER_EBT || method == SOLVER_IEBT || method == SOLVER_ABM){ // EBT and ABM have very simular structure so use the same density calculation algo
		double xm = spp->getX(0)[dim]+1e-6;

		std::vector<point> points(breaks.size());

		// int current_interval = breaks.size()-2;
		for (int i=0; i<spp->J; ++i){ // loop over all cohorts 
			double x = spp->getX(i)[dim];
			double N = spp->getU(i);
			if (method == SOLVER_EBT || method == SOLVER_IEBT) if (i == spp->J-1) x = spp->xb[dim] + x/(N+1e-12); // For EBT, real-ize x if it's boundary cohort
			
			int current_interval = std::upper_bound(breaks.begin(), breaks.end(), x) - breaks.begin() -1; // this assumes breaks is sorted ascending
			current_interval = std::clamp<int>(current_interval, 0, breaks.size()-1);
			// while(breaks[current_interval]>x) --current_interval; // decrement interval until x fits
			// cout << current_interval << ", x = " << x << "(" << N << ") in [" << breaks[current_interval] << ", " << breaks[current_interval+1] << "]\n"; cout.flush();
			if (N>0){
				points[current_interval].abund += N;
				points[current_interval].count += 1;
				points[current_interval].xmean += N*x;
			}
		}
		// Compute mean x in each interval (each point corresponds to 1 interval)
		for (int i=0; i<points.size(); ++i) if (points[i].count>0) points[i].xmean /= points[i].abund;
		
		// remove 0-count points (i.e., delete intervals with no cohorts)
		auto pred = [this](const point& p) -> bool {
			return p.count == 0; 
		};

		auto pend = std::remove_if(points.begin(), points.end(), pred);
		points.erase(pend, points.end());

		// cout << "mean x and abund (removed):\n";
		// for (int i=0; i<points.size(); ++i) cout << i << "\t" << points[i].count << "\t" << points[i].xmean << "\t" << points[i].abund << "\n";	
		// cout << "--\n";

		// Now treat xmean as the x std::vector and calculate the width spanned by each point, i.e., 1D Voronoi cells
		// to get u, divide abundance by width for each point
		if (points.size() > 2){
			std::vector<double> h(points.size());
			h[0] = (points[1].xmean+points[0].xmean)/2 - spp->xb[dim];
			for (int i=1; i<h.size()-1; ++i) h[i] = (points[i+1].xmean - points[i-1].xmean)/2;
			h[h.size()-1] = xm - (points[h.size()-1].xmean+points[h.size()-2].xmean)/2;

			xx.reserve(points.size());
			uu.reserve(points.size());
			for (int i=0; i<points.size(); ++i){
				xx.push_back(points[i].xmean);
				uu.push_back(points[i].abund / h[i]);
			}
		}

	}

	// Note: This only works for 1D models 
	else if (method == SOLVER_CM || method == SOLVER_ICM){
		xx.reserve(spp->J);
		uu.reserve(spp->J);
		for (int i=spp->J-1; i>=0; --i){
			xx.push_back(spp->getX(i)[dim]);
			uu.push_back(spp->getU(i));
		}
		
	}	
	
	// in grid-based methods, we aggregate all dimensions other then dim
	else {
		std::vector<point> points(spp->X[dim].size());

		// int current_interval = breaks.size()-2;
		for (int i=0; i<spp->J; ++i){ // loop over all cohorts 
			double x = spp->getX(i)[dim];
			double U = spp->getU(i);
			
			std::vector<double> dx = utils::tensor::coord_value(utils::tensor::index(i, spp->dim_centres), spp->h);
			double dV = std::accumulate(dx.begin(), dx.end(), 1.0, std::multiplies<double>()); 

			int current_interval = std::upper_bound(spp->X[dim].begin(), spp->X[dim].end(), x) - spp->X[dim].begin() -1; // this assumes breaks is sorted ascending
			current_interval = std::clamp<int>(current_interval, 0, points.size()-1);

			// cout << current_interval << ", x = " << x << "(" << N << ") in [" << breaks[current_interval] << ", " << breaks[current_interval+1] << "]\n"; cout.flush();
			points[current_interval].abund += U*dV;
			points[current_interval].count += 1;
			// points[current_interval].h     += spp->h[dim][current_interval];
			points[current_interval].xmean += U*dV*x;
		}

		if (points.size() > 2){
			xx.reserve(points.size());
			uu.reserve(points.size());
			for (int i=0; i<points.size(); ++i){
				xx.push_back(points[i].xmean/points[i].abund);
				uu.push_back(points[i].abund / spp->h[dim][i]); // /(points[i].h/points[i].count));
			}
		}

	}	


	// finally, interpolate to values in breaks
	if (xx.size() >= 2){ 
		
		Spline spl;
		spl.splineType = Spline::LINEAR; //Spline::CONSTRAINED_CUBIC;
		spl.extrapolate = extrapolation_method; //Spline::ZERO; //Spline::ZERO;
		spl.set_points(xx, uu);
		 
		std::vector <double> dens;
		dens.reserve(xx.size());
		for (int i=0; i<breaks.size(); ++i){
			dens.push_back(spl.eval(breaks[i]));			
		}
		
		return dens;
	}
	else {
		return std::vector<double>(breaks.size(), 0);
	}
}




// // breaks must be sorted ascending
// std::vector<std::vector<double>> Solver::getDensitySpecies2D(int k, const std::vector<int>& axes, const std::vector<std::vector<double>>& breaks, Spline::Extr extrapolation_method){
// 	auto spp = species_vec[k];
// 	std::vector <double> xx, uu;

// 	if (axes.size() != 2) throw std::runtime_error("axes should contain 2 elements.");
// 	if (breaks.size() != axes.size()) throw std::runtime_error("'breaks' should be a vector of n vectors, where n is the number of axes specified via the 'axes' argument.");
	
// 	std::vector<int> dims = {(int)breaks[0].size(), (int)breaks[1].size()};

// 	// if (method == SOLVER_ABM) spp->sortCohortsDescending(dim); // ABM needs sorted cohorts
// 	// if (method == SOLVER_EBT || method == SOLVER_IEBT) spp->sortCohortsDescending(dim, 1); // for EBT, skip boundary cohort and sort

// 	if (method == SOLVER_EBT || method == SOLVER_IEBT || method == SOLVER_ABM){ // EBT and ABM have very simular structure so use the same density calculation algo

// 		realizeEbtBoundaryCohort(spp);
// 		std::vector<point> points(breaks.size());

// 		// int current_interval = breaks.size()-2;
// 		for (int i=0; i<spp->J; ++i){ // loop over all cohorts 
// 			std::vector<double> x = spp->getX(i);
// 			double N = spp->getU(i);

// 			std::vector<int> current_interval(axes.size()); 
// 			for (int j=0; j<axes.size(); ++j){
// 				current_interval[j] = std::upper_bound(breaks[j].begin(), breaks[j].end(), x[axes[j]]) - breaks[j].begin() -1; // this assumes breaks is sorted ascending
// 				current_interval[j] = std::clamp<int>(current_interval[j], 0, breaks[j].size()-1);
// 			}
// 			// while(breaks[current_interval]>x) --current_interval; // decrement interval until x fits
// 			// cout << current_interval << ", x = " << x << "(" << N << ") in [" << breaks[current_interval] << ", " << breaks[current_interval+1] << "]\n"; cout.flush();
// 			if (N>0){
// 				int loc = utils::tensor::location(current_interval, dims);
// 				points[loc].abund += N;
// 				points[loc].count += 1;
// 				points[loc].xmean_vec.resize(axes.size());
// 				for (int j=0; j<axes.size(); ++j){
// 					points[loc].xmean_vec[j] += N*x[axes[j]];
// 				}
// 			}
// 		}
// 		// Compute mean x in each interval (each point corresponds to 1 interval)
// 		for (int i=0; i<points.size(); ++i){
// 			if (points[i].count>0){
// 				for (int j=0; j<axes.size(); ++j){
// 					points[i].xmean_vec[j] /= points[i].abund;
// 				}
// 			}
// 		}
// 		// remove 0-count points (i.e., delete intervals with no cohorts)
// 		auto pred = [this](const point& p) -> bool {
// 			return p.count == 0; 
// 		};

// 		auto pend = std::remove_if(points.begin(), points.end(), pred);
// 		points.erase(pend, points.end());

// 		// cout << "mean x and abund (removed):\n";
// 		// for (int i=0; i<points.size(); ++i) cout << i << "\t" << points[i].count << "\t" << points[i].xmean << "\t" << points[i].abund << "\n";	
// 		// cout << "--\n";

// 		restoreEbtBoundaryCohort(spp);

// 	}

// 	else if (method == SOLVER_CM || method == SOLVER_ICM){
// 		xx.reserve(spp->J);
// 		uu.reserve(spp->J);
// 		for (int i=spp->J-1; i>=0; --i){
// 			xx.push_back(spp->getX(i)[dim]);
// 			uu.push_back(spp->getU(i));
// 		}
		
// 	}	
	
// 	else {
// 		// in grid-based methods, we must aggregate all dimensions other then dim
// 		xx.reserve(spp->J);
// 		uu.reserve(spp->J);
// 		for (int i=0; i<spp->J-1; ++i){
// 			xx.push_back(spp->getX(i)[dim]);
// 			uu.push_back(spp->getU(i));
// 		}
		
// 	}	


// 	// interpolate density at each value in breaks from uu
// 	if (xx.size() >= 2){ 
		
// 		Spline spl;
// 		spl.splineType = Spline::LINEAR; //Spline::CONSTRAINED_CUBIC;
// 		spl.extrapolate = extrapolation_method; //Spline::ZERO; //Spline::ZERO;
// 		spl.set_points(xx, uu);
		 
// 		std::vector <double> dens;
// 		dens.reserve(xx.size());
// 		for (int i=0; i<breaks.size(); ++i){
// 			dens.push_back(spl.eval(breaks[i]));			
// 		}
		
// 		return dens;
// 	}
// 	else {
// 		return std::vector<double>(breaks.size(), 0);
// 	}
// }
