#ifndef PSPM_ODE_RKCK45_H_
#define PSPM_ODE_RKCK45_H_

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 *  Code for RK Solvers 
 *  Original version obtained from:
 *     https://www.physics.rutgers.edu/grad/509/src_numerics/ODE/1/ode1.cc 
 *  
 *  v2
 *  Modified sligtly by Jaideep Joshi (30/5/2020)
 *  - Cleaned up indentation
 *  - Better class names
 *  - remove requirement of stop time
 *  - Add "step_to" function
 *  - move all containers to main class for ease of resizing
 *
 *  v3
 *  - uses GSL style step change algorithm
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <cassert>
#include <fstream>


//using namespace std;

typedef std::vector<double> container;

//template <class functor, class container>
//void Euler(double x, double h, container& y, functor& derivs){
	//container fk(y.size());
	//derivs(x, y, fk);
	//for (int i=0; i<y.size(); i++) y[i] += h*fk[i]; 
//}

/////////////////////////////////////////////////////////////////////////////
///////////////// METHOD OF RUNKGE-KUTTA 4-th ORDER   ///////////////////////
/////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
///////// METHOD OF RUNGE-KUTTA 5-th ORDER, ADAPTIVE STEP  /////////////
////////////////////////////////////////////////////////////////////////
class RKCK45{
	private:  
	static constexpr double SAFETY = 0.9;
	static constexpr double PGROW  = -0.2;
	static constexpr double PSHRNK = -0.25;
	static constexpr double ERRCON = 1.89e-4;
	static constexpr double a_y = 1;
	static constexpr double a_dydt = 0;
	container yscal;   // scaling factors
	double ht;         // current time-step
	double eps_rel, eps_abs;        // Accuracy we check at each step
	double xt, t_stop; // start and stop time
	int nok=0, nbad=0; // number of bad and good steps
	container dydx;    // temporary current derivatives
	int nfe = 0;
	double hmin = 1e-6;

	container k1, k2, k3, k4, k5, yt;
	
	int sys_size = 0;

	public:
	// Constructor takes the following arguments:
	// 1) t_start  -- start time          
	// 2) t_stop   -- end_time            
	// 3) size     -- number of equations 
	// 4) accuracy -- desired accuracy
	// 5) h1       -- trial size of the first step
	RKCK45(double t_start_, double accuracy, double h1);
	RKCK45(double t_start_, double accuracy, double h1, double _hmin);
	~RKCK45();

	// Resize the container 
	void resize(int new_size);
	
	// The main function which takes next RK step with predefined accuracy
	// returns bool:   true if time becomes larger of equal to t_stop, otherwise false
	// arguments:
	//           1) x      -- current time (independent variable), which will be changes for a step when successful
	//           2) y      -- current solution (dependent variable)
	//           3) derivs -- function which evaluates derivatives
	template <class functor>
	void Step(double& x, container& y, functor& derivs, double hmax=1e20);
 
	//template <class functor>
	//void Step_RK4(double &x, double h, container& y, functor& derivs);
	
	template <class functor, class AfterStep>
	void Step_to(double t_stop, double& x, container& y, functor& derivs, AfterStep &after_step);

	double h(){return ht;}
	double size(){return sys_size;}
	int get_fn_evals(){return nfe;}

	void save(std::ostream &fout);
	void restore(std::istream &fin);
	
	private:
	template <class functor>
	void RKTry(container& y, container& dydx, double& x, double h, container& yout, container& yerr, functor& derivs);
	
	template <class functor>
	void RKStep(container& y, container& dydx, double& x, double htry, double& hdid, double& hnext, functor& derivs);

};


#include "../src/rkck45.tpp"


#endif



