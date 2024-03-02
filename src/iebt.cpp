#include <cassert>
#include <string>
#include <algorithm>
// #include <execution>
#include "solver.h"
#include "linear_system.h"
using namespace std;

// void Solver::calcRates_iEBT(double t, vector<double>::iterator S, vector<double>::iterator dSdt){
// 	vector<double>::iterator its = S    + n_statevars_system; // Skip system variables
// 	vector<double>::iterator itr = dSdt + n_statevars_system;
	
// 	for (int s = 0; s<species_vec.size(); ++s){
// 		auto spp = species_vec[s];
		
// 		its += (n_statevars_internal)*spp->J; // skip x and u 
// 		for (int i=0; i<n_statevars_internal*spp->J; ++i) *itr++ = 0; // set dx/dt and du/dt to 0 
	
// 		if (spp->n_extra_statevars > 0){
// 			auto itr_prev = itr;
// 			spp->getExtraRates(itr); // TODO/FIXME: Does calc of extra rates need t and env?
// 			assert(distance(itr_prev, itr) == spp->n_extra_statevars*spp->J);
// 			its += spp->n_extra_statevars*spp->J; 	
// 		}
// 	}

// }


double identity_matrix(int i, int j){
	return (i==j)? 1:0;
}

// Note: Boundary cohort dynamics for IEBT can be written in matrix form:
//  /            [  0   -dmu/dx0      ... -dmu/dxn ]  \    [ N0(t+1) ] = [ N0 + B*dt ]      
// |             [  g0   dg0/dx0      ...  dg0/dxn ]   |   [ p0(t+1) ] = [ p0        ]   
// |(1+mu*dt)I - [  g1   dg1/dx0      ...  dg1/dxn ]dt | * [ p1(t+1) ] = [ p1        ]               
// |             [  ..   ..           ...   ..     ]   |   [ ..      ] = [ ..        ]   
//  \            [  gn   dgn/dx0      ...  dgn/dxn ]  /    [ pn(t+1) ] = [ pn        ]   
//               <---------------- A -------------->                     <---- B ---->
//  /                         \.     
// |(1+mu*dt)I - [ 0   -mx ]dt | X(t+1) = X(t)+ Bdt*[1]
// |             [ g'   gx']   |                    [0]
//  \                         /
void Solver::stepU_iEBT(double t, vector<double> &S, vector<double> &dSdt, double dt){
		
	vector<double>::iterator its = S.begin() + n_statevars_system; // Skip system variables
	
	// 1. Take implicit step for U
	for (int s = 0; s<species_vec.size(); ++s){
		Species_Base* spp = species_vec[s];
		
		// [S S S x u x u x u x u x u a b c a b c a b c a b c a b c] <--- full SV for species is this
		//        ^ its, itr
		double *XU = &(*its); 
		int J = spp->J;

		std::vector<double> pi0 = spp->getX(spp->J-1);	 // last cohort is pi0, N0
		double N0 = spp->getU(spp->J-1);
		// std::cout << "Species " << s << ", t = " << t << ", pi = " << pi0 << ", N0 = " << N0 << "\n";

		std::vector<double> grad_dx(spp->istate_size, control.ebt_grad_dx);

		std::vector<std::vector<double>> g_gx = spp->growthRateGradient(-1, t, env, grad_dx);
		std::vector<double> m_mx = spp->mortalityRateGradient(-1, t, env, grad_dx);
		// std::cout << "Species " << s << ", t = " << t << ", ";
		// std::cout << "g = " << g_gx[0] << ", gx = " << g_gx[1] << "\n";
		// std::cout << "Species " << s << ", t = " << t << ", ";
		// std::cout << "m = " << m_mx[0] << ", mx = " << m_mx[1] << "\n";

		double birthFlux;
		double pe = spp->establishmentProbability(t, env);
		if (spp->birth_flux_in < 0){	
			birthFlux = calcSpeciesBirthFlux(s,t) * pe;
		}
		else{
			birthFlux = spp->birth_flux_in * pe;
		}
		if (debug) std::cout << "Species " << s << ", t = " << t << ", bf = " << birthFlux << '\n';


		// internal cohorts
		// for (int i=0; i<J-1; ++i){
		// 	XU[2*i+0] += spp->growthRate(i, spp->getX(i), t, env)*dt;
		// 	XU[2*i+1] /= 1+spp->mortalityRate(i, spp->getX(i), t, env)*dt;
		// }
		// below version allows parallel execution
		vector<int>a(J-1);
		std::iota(a.begin(), a.end(), 0);
		std::for_each(
//			std::execution::par,
			a.begin(),
			a.end(),
			[this, XU, spp, t, dt](int i){
				int nk = n_statevars_cohort(spp);
				std::vector<double> g = spp->growthRate(i, t, env);	 // dx/dt
				for (int k=0; k<spp->istate_size; ++k){
					XU[nk*i+k] += g[k]*dt;
				}
				// std::cout << "   spp = " << spp << ", x = " << spp->getX(i) << " , mort = " << spp->mortalityRate(i, t, env) << "\n";
				XU[nk*i+spp->istate_size] /= 1+spp->mortalityRate(i, t, env)*dt;
			}
		);

		// pi0 cohort
		int K = spp->istate_size;
		Matrix A(K+1, Vector(K+1, 0)); // create (n+1)x(n+1) Zero matrix

		// 1. Construct the matrix -A*dt (note negative A)
		// top-left element is 0
		A[0][0] = 0; // included for clarity. 
		// first row (i=0), columns j=[1..K] contains mortality gradient vector
		for (int j=0; j<K; ++j){
			A[0][j+1] = m_mx[j+1]*dt; // m_mx[1..K] is mortality gradient, m_mx[0] is the mortality rate
		}
		// all other rows (i=1..K), columns (j=0..K) contain the transposed g_gx matrix
		for (int i=0; i<K; ++i){ // row index goes from 0..K-1
			for (int j=0; j<K+1; ++j){ // column index goes from 0..K
				A[i+1][j] = -g_gx[j][i]*dt; 
				// ^ i+1 here because we need to skip the first row in A
			}
		}
		// 2. Add (1+mu*dt)I to (-A*dt)
		for (int i=0; i<K+1; ++i){
			A[i][i] += 1 + m_mx[0]*dt;
		}

		// 3. Construct B
		Vector B(K+1, 0);
		B[0] = N0 + birthFlux*dt;
		for (int i=0; i<K; ++i){
			B[i+1] = pi0[i];
		}
		// std::cout << "Species " << s << ", t = " << t << ", a1/b1/c1/a2/b2/c2 = " << A[0][0] << "/" << A[0][1] << "/" << B[0] << "/" << A[1][0] << "/" << A[1][1] << "/" << B[1] << '\n';

		Vector X = lupSolve(A, B);

		// double a1 = 1 + mb*dt;
		// double b1 = mortGrad*dt;
		// double c1 = birthFlux*dt + N0;
		// double a2 = -gb*dt;
		// double b2 = 1 - growthGrad*dt + mb*dt;
		// double c2 = pi0;

		// double pinew = (a2*c1-a1*c2)/(a2*b1-a1*b2);
		// double unew = (b2*c1-b1*c2)/(b2*a1-b1*a2);
		
		int nk = n_statevars_cohort(spp);
		for (int k=0; k<spp->istate_size; ++k){
			if (X[k+1] < 0) throw std::runtime_error("pi_"+std::to_string(k)+" < 0: "+std::to_string(X[k+1]));
			XU[nk*(J-1)+k] = X[k+1];
		}
		if (X[0] < 0) throw std::runtime_error("u0 < 0: "+std::to_string(X[0]));
		XU[nk*(J-1)+spp->istate_size] = X[0];
		// std::cout << "Species " << s << ", t = " << t << ", N0 = " << X[0] << ", pi0 = " << X[1] << '\n';
		
		its += J*(nk+spp->n_accumulators);

	}
	
	assert(distance(S.begin(), its)==S.size());
}


