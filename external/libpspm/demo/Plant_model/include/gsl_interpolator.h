// -*-c++-*-
#ifndef PLANT_PLANT_GSL_INTERPOLATOR_H_
#define PLANT_PLANT_GSL_INTERPOLATOR_H_

#include <iostream>
#include <cassert>
#include <cmath>
#include <list>
#include <string>
#include <vector>
#include <fstream>
#include <iomanip>
#include <chrono>

#include <gsl/gsl_spline.h>
//using namespace std::chrono;


class AdaptiveInterpolator{
	private:
	gsl_spline * spline = NULL;
	gsl_interp_accel  * accel = NULL;

	int npoints_0 = 17;
	int max_depth = 16;
	
	std::string type = "cspline";

	int npoints; 
	double dx, dxmin;
	
	double rel_tol = 1e-4;
	double abs_tol = 1e-4;
	double a, b;
	
	public:
	std::list<double> xx, yy; 
	std::list<bool> zz;
	
	public:
	AdaptiveInterpolator(){
		npoints = npoints_0;
		accel = gsl_interp_accel_alloc();
	}
	
	
	~AdaptiveInterpolator(){
		gsl_spline_free(spline);
		gsl_interp_accel_free(accel);
	}


	inline void realloc(){
		gsl_spline_free(spline);
		if      (type == "cspline") spline = gsl_spline_alloc(gsl_interp_cspline, npoints);
		else if (type == "akima")   spline = gsl_spline_alloc(gsl_interp_akima,   npoints);
		else { assert(false && "Invalid interpolation method"); }
		gsl_interp_accel_reset(accel);
	}


	template <class T>
	void constructAndReset(std::vector <T> &X, std::vector<T> &Y){
		xx = std::list<T>(X.begin(), X.end());
		yy = std::list<T>(Y.begin(), Y.end());
		construct(X,Y);
	}
	

	template <class T>
	void construct(std::vector <T> &X, std::vector<T> &Y){
		assert(X.size() == Y.size());
		npoints = X.size();
		
		assert(npoints > 2);
					
		realloc();
		gsl_spline_init(spline, X.data(), Y.data(), X.size());
	}

		
	template <typename Function>  
	void construct(Function f, double a, double b){

//		auto t1 = std::chrono::steady_clock::now();
		assert(a < b);
		
		realloc();	
				
		npoints = npoints_0;		
		dx = (b - a)/(npoints - 1);
		dxmin = dx / pow(2, max_depth);
		
		xx.clear(); yy.clear(); zz.clear();
		for (int i = 0; i < npoints; ++i) {
			const double xi = a + i*dx;
			xx.push_back(xi);
			yy.push_back(f(xi));
			zz.push_back(true);
		}
		
		std::vector<double> X(xx.begin(), xx.end());
		std::vector<double> Y(yy.begin(), yy.end());
		// set gsl spline here
		construct(X, Y);		

//		using namespace std;
//		cout << "Constructed light profile:\n";
//		for (int i=0; i<X.size(); ++i) cout << X[i] << " " << Y[i] << "\n";
//		cout << "~~~~~~~~~~~~~~~" << endl;		

		bool b_error = true;
		while (b_error) {
			b_error = refine(f);
		}

//		auto t2 = std::chrono::steady_clock::now();
//		std::cout << "Interpolator Construction: (" << xx.size() << " points) [" << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() << " usec]" << std::endl;

	}
	

	inline double eval(double x) const{
//		std::cout << std::setprecision(12) <<  "Interpolator::eval( " << xx.front() << " --- " << x << " --- " << xx.back() << std::endl; 
		if (x <= xx.front()) return yy.front();
		if (x >= xx.back())  return yy.back();
 		return gsl_spline_eval(spline, x, accel);
	}


	inline bool check_no_err(double y_true, double y_pred) const {
		const double err_abs = fabs(y_true - y_pred);
		const double err_rel = fabs(1 - y_pred / y_true);
		return err_abs < abs_tol || err_rel < rel_tol;	// FIXME: abs_tol doesnt make sense! What if you're interpolating {1e-22, 2e-22, 1e-22} to {1e-22, 2e-8, 1e-22}? -- NO, it does! 
	}
	

//	template <typename Function>
//	bool refine(Function f) {
//	  dx /= 2;

//	  assert(dx >= dxmin && "Message: Maximum refinement reached");

//	  bool flag = false;

//	  std::list<double>::iterator xi = ++xx.begin(), yi = ++yy.begin();
//	  std::list<bool>::iterator zi = ++zz.begin();
//	  for (; xi != xx.end(); ++xi, ++yi, ++zi) {
//		if (*zi) {
//		  const double x_mid = *xi - dx;
//		  const double y_mid = f(x_mid);
//		  const double p_mid = eval(x_mid);

////			std::cout << "xm,ym,pm = " << x_mid << ", " << y_mid << ", " << p_mid << std::endl; 

//		  // Flag for refinement based on
//		  const bool mid_error = !check_no_err(y_mid, p_mid);
//		  *zi = mid_error;
//		  if (mid_error){	// FIXME: Need to rigorously test this logic
//			  xx.insert(xi, x_mid);
//			  yy.insert(yi, y_mid);
//			  zz.insert(zi, mid_error);
//			  ++npoints;
//		  }
//			
//		  flag = flag || mid_error;
//		}
//	  }
//	  
////	  std::cout << "npoints = " << npoints << std::endl;
//	  
//		std::vector<double> X(xx.begin(), xx.end());
//		std::vector<double> Y(yy.begin(), yy.end());
//		construct(X, Y);		
//	  
//		return flag;
//	}
	
// *** ORIGINAL REFINING CODE FROM PLANT, with iterators ++ed to prevent extrapolation
	template <typename Function>
	bool refine(Function target) {
	  dx /= 2;

	  if (dx < dxmin) {
		assert(false);
	  }

	  bool flag = false;

	  std::list<double>::iterator xi = ++xx.begin(), yi = ++yy.begin();
	  std::list<bool>::iterator zi = ++zz.begin();
	  for (; xi != xx.end(); ++xi, ++yi, ++zi) {
		if (*zi) {
		  const double x_mid = *xi - dx;
		  const double y_mid = target(x_mid);
		  const double p_mid = eval(x_mid);

//			std::cout << "xm,ym,pm = " << x_mid << ", " << y_mid << ", " << p_mid << std::endl; 

		  // Always insert the new points.
		  xx.insert(xi, x_mid);
		  yy.insert(yi, y_mid);

		  // Flag for refinement based on
		  const bool flag_mid = !check_no_err(y_mid, p_mid);
		  // If error was OK/not OK (flag_mid is true/false), say that
		  // this interval is OK/not OK...
		  *zi = flag_mid;
		  // ...and that that the interval implied by the mid point is
		  // OK/not OK.
		  zz.insert(zi, flag_mid);

		  flag = flag || flag_mid;
		}
	  }

//	  std::cout << "Refined with npoints = " << npoints << std::endl;

	  // Recompute interpolator to use new points added during refinement.
		std::vector<double> X(xx.begin(), xx.end());
		std::vector<double> Y(yy.begin(), yy.end());
		construct(X, Y);		

	  return flag;
	}
	
	
	void print(std::ostream & fout = std::cout){
//		std::ofstream fout(filename);
		fout << std::setprecision(12);
		for (auto x : xx) fout << x << "\t";
		fout << std::endl;
		for (auto y : yy) fout << y << "\t";
		fout << std::endl;
		fout.flush();
//		fout.close();
	}
	
	
};







#endif

