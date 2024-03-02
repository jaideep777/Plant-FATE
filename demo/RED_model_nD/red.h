#ifndef DEMO_RED_MODEL_H
#define DEMO_RED_MODEL_H

#include <individual_base.h>

class LightEnvironment : public EnvironmentBase{
	
	double E = 0;

	public:
	double evalEnv(double xn, double t){
		return E;
	}

	// This function must do any necessary precomputations to facilitate evalEnv()
	// Therefore, this should calculate env for all X when it is a function of X
	// In such a case, the solver's SubdivisionSpline can be ussed
	// Note: The state vector in the solver will not be updated until the RK step is completed. 
	// Hence, explicitly pass the state to this function.
	void computeEnv(double t, Solver * S, std::vector<double>::iterator s, std::vector<double>::iterator dsdt){
		//             _xm 
		// Calculate _/ w(z,t)u(z,t)dz
		//         xb
		auto w = [S](int i, double t) -> double {
			std::vector<double> z = S->species_vec[0]->getX(i);
			return 0.396*pow(z[0], 0.749)/10000;
		};
		E = S->state_integral(w, t, 0);
	}

};



class RED_Plant : public IndividualBase<2>{
	public:

	//double input_seed_rain = 1;	
	vector <string> varnames;

	int nrc = 0; // number of evals of compute_vars_phys() - derivative computations actually done by plant
	int ndc = 0; // number of evals of mortality_rate() - derivative computations requested by solver
	int nbc = 0;

	double a0 = 0.396;
	double phiA = 0.749;
	double g0 = 0.0838;
	double phiG = 0.7134;
	double m0 = 1;
	double mort = 0.035;
	double alpha = 0.1;
	double beta = 0.002;
	double mu0;

	RED_Plant() {
		mu0 = mort*m0/g0;
	}

	void set_size(const std::array <double, 2>& _x){
	}

	double init_density(void * env, double input_seed_rain){
		return 100/pow(x[0],4);
	}

	void preCompute(double t, void * env){
	}

	double establishmentProbability(double t, void * env){
		return 1;
	}

	std::array<double,2> growthRate(double t, void * env){
		++nrc;
		// std::cout << "growth Rate Gradient compute" << std::endl;
		// std::cout << xn << std::endl;
		double x0 = g0*pow(x[0],phiG);
		// std::cout << x0 << std::endl;	
		double x1 = -beta * x[1];
		// std::cout << x1 << std::endl;
		// std::cout <<x0 << "   " <<  x1 << std::endl;
		// std::cout << "growth Rate Gradient compute: END" << std::endl;

		return {x0, x1};
	}

	double mortalityRate(double t, void * env){
		++ndc;
		return mort;
	}
	
	double birthRate(double t, void * env){
		++nbc;
		if (x[0] < 0) throw std::runtime_error("x_0 is negative");
		LightEnvironment* env1 = (LightEnvironment*)env;
		// std::cout << "Birth Rate for " << x[0] << std::endl;
		return 0.1/0.9*g0*pow(x[0],phiG)*(1-env1->evalEnv(x[0],t));
	}

	void init_accumulators(double t, void * env){
	}

	vector<double>::iterator set_accumulators(vector<double>::iterator &it){
		return it;
	}
	vector<double>::iterator get_accumulators(vector<double>::iterator &it){
		return it;
	}
	vector<double>::iterator get_accumulatorRates(vector<double>::iterator &it){
		return it;
	}

	void print(std::ostream& out = std::cout){
	}
};


#endif
