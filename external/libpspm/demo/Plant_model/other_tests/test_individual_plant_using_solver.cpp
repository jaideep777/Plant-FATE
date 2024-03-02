#include <vector>
#include <fstream>
#include "../include/plant.h"
#include "../include/gsl_ode_solver.h"

//#pragma GCC diagnostic push 
//#pragma GCC diagnostic ignored "-Wterminate"
//#include "../third_party_includes/gnuplot_i.hpp"
//#pragma GCC diagnostic pop

using namespace std;
using namespace plant;

// This code reproduces the growth of individual plants (1st Plot) from https://traitecoevo.github.io/plant/articles/individuals.html

int main(){
	
	initPlantParameters(par);

	Environment env(1);
//	env.canopy_openness = fixed_canopy_openness;
	cout << "Env: " << env.canopy_openness(1) << endl;
	
	
	Plant p;
	p.lma = 0.1978791;
	p.set_height(0.3441948);
	par.r_l   = 198.4545; //39.27 / 0.1978791; // JAI: Should be 39.27/lma; 
	par.k_l = 0.4565855;
	
	cout << p << endl;
	
	p.compute_vars_phys(env, false);

	cout << p << endl;
	
	
	auto func_lambda = [&p, &env](double t, const double y[], double f[]) -> int {
		env.time=t;
		p.set_state(y);
		p.compute_vars_phys(env, false);
		p.get_rates(f);
		
		return GSL_SUCCESS;
	};
	
	
	ODE_Solver<decltype(func_lambda)> sys(func_lambda, 5);

	vector <double> state0 = {p.vars.height, p.vars.mortality, p.vars.fecundity, p.vars.area_heartwood, p.vars.mass_heartwood};

	sys.initialize(state0);

	double t0 = 0, tf = 50, dt = 1;
	size_t nsteps = (tf-t0)/dt;
//	vector <double> heights;// = {p.vars.height};

	ofstream fout("ind_plant.txt");

	for (size_t i=0; i < nsteps; ++i){

		sys.step_to(i*dt);		
		
		fout << i*dt << "\t" <<
			p.vars.height         << "\t" <<
			p.vars.mortality      << "\t" <<
			p.vars.fecundity      << "\t" <<
			p.vars.area_heartwood << "\t" <<
			p.vars.mass_heartwood << endl;
						
//		cout << p.vars.fecundity << " ";		
//		heights.push_back(p.vars.height);
	}
	
	cout << p << endl;

//	Plant p2;
//	p2.step_by(50, env);
//	cout << p2 << endl;
	
//	for (size_t i=0; i<heights.size(); ++i) cout << heights[i] << " ";
//	cout << endl;	
	
//	Gnuplot g1("lines");
//    g1.set_style("impulses").plot_x(heights,"user-defined doubles");


//	for (auto i : heights) fout << i << endl;

//	int n; cin >> n;
	
	return 0;
}


