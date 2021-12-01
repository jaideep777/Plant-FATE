#ifndef MATH_UTILS_LAMBERT_W_H_
#define MATH_UTILS_LAMBERT_W_H_

#include <iostream>
#include <cassert>
#include <cmath>

inline double lambertw0(double x, double tol=1e-8, int maxiter = 100){
	assert(x>=0);
	if (x == 0) return 0;

	int niter = 0;
	double err = 1;
	double w = (x < 2.718282)? 0.4*x : log(x) - log(log(x)) + 0.3;  // first approximation - works great for 0 < x < 100 - causes convergence in 2 iterations

	while(niter < maxiter){
		double ew = exp(w);
		double wew = w*ew;
		double r = wew-x; // Residual, Eq. 5.2
		if (fabs(r) < tol) break;
		
		// Halley's method, for Lambert W from Corless, et al. 1996, Eq. 5.9
		w = w-r/(wew+ew-0.5*r*(w+2)/(w+1));
		++niter;
	}

	if (niter >= maxiter) std::cout << "Warning: lambertw0() - max iterations reached. returning last approximation\n";	
	std::cout << "LambertW: x = " << x << ", w(x) = " << w << ",iter = " << niter << "\n";
	return w;
}

#endif

