#ifndef __PSPM_TEST_TEST_TEST_MODEL_H_
#define __PSPM_TEST_TEST_TEST_MODEL_H_

#include <individual_base.h>
#include <environment_base.h>
#include <solver.h>

class Environment : public EnvironmentBase{
	
	double E = 0;

	public:
	double evalEnv(double x, double t){
		return E;
	}

	void computeEnv(double t, Solver * S, std::vector<double>::iterator s, std::vector<double>::iterator dsdt){
		//             _xm 
		// Calculate _/ w(z,t)u(z,t)dz
		//         xb
		auto w = [S](int i, double t) -> double {
			double z = S->species_vec[0]->getX(i)[0];
			if (z <= 1.0/3) 
				return 1;
			else if (z > 1.0/3 && z <= 2.0/3) 
				return pow(2-3*z, 3)*(54*z*z-27*z+4);
			else 
				return 0;
		};
		E = S->state_integral(w, t, 0);
	}

};



class Plant{
	public:
	double height;
	double mortality = 0;
	double viable_seeds = 0;
	double heart_mass = 0;
	double sap_mass = 0;

	
	Plant(double h){
		height = h;
	}

	std::vector<double> calcRates(){
		std::vector<double> rates(4);
		rates[0] = -2;
		rates[1] = 0;
		rates[2] = -20*height;
		rates[3] = -30*height;
		return rates;
	}
};





class TestModel : public Plant, public IndividualBase<1>{
	public:
	double sc = 10;

	// precomputed demographic rates
	double g,m,f,se;
	
	TestModel() : Plant(0) {}

	double init_density(void * _env, double bf){
		double _x = x[0];
		return pow(1-_x,2)/pow(1+_x,4) + (1-_x)/pow(1+_x,3);
	}
	
	void preCompute(double t, void * _env){
		Environment* env = (Environment*)_env;
		double _x = x[0];	
		// growth rate
		double E = env->evalEnv(_x,t);
		double a = 0.16+0.22*exp(-0.225*t*t);
		g = {0.225*(1-_x*_x)*(E/(1+E*E))*t*(1+a*a)/a};

		m = 1.35*t*E/a;
		
		double oneplusa = 1+a;
		double n1 = 0.225*t*_x*_x*(1-_x)*(1-_x)*E/(1+E)/(1+E)*oneplusa*oneplusa/a;
		double n2 = (1+exp(-0.225*t*t))/(61-88*log(2)+(38*log(2)-79.0/3)*exp(-0.225*t*t));
		f = n1*n2;

		se = 1;

	}

	std::array<double,1> growthRate(double t, void * _env){
		return {g};
	}

	double mortalityRate(double t, void * _env){
		return m;
	}

	double birthRate(double t, void * _env){
		return f;
	}

	double establishmentProbability(double t, void  * _env){
		return se;
	}


	void set_size(const std::array <double, 1>& _x){
		height = _x[0];
		mortality = 0.1*height + 1e-12; //exp(-height);	
		viable_seeds = 100*height + 1e-13;
		heart_mass = 1000*height + 1e-14;
		sap_mass = 10*height + 1e-15;
	}

	void init_accumulators(double t, void * env){
	}

	std::vector<double>::iterator set_accumulators(std::vector<double>::iterator &it){
		mortality    = *it++;
		viable_seeds = *it++;
		heart_mass   = *it++;
		sap_mass     = *it++;
		return it;
	}

	std::vector<double>::iterator get_accumulators(std::vector<double>::iterator &it){
		*it++ = mortality;
		*it++ = viable_seeds;
		*it++ = heart_mass;
		*it++ = sap_mass;
		return it;
	}

	std::vector<double>::iterator get_accumulatorRates(std::vector<double>::iterator &it){
		std::vector<double> r = calcRates();
		*it++ = r[0];
		*it++ = r[1];
		*it++ = r[2];
		*it++ = r[3];
		return it;
	}

	void print(std::ostream &out = std::cout){
		out << std::setw(10) << mortality << "\t"
		    << std::setw(10) << viable_seeds << "\t"
			<< std::setw(10) << heart_mass << "\t"
			<< std::setw(10) << sap_mass << "\t";
	}

	void save(std::ostream& fout){
		fout << "TestModel::v1 ";
	}

	void restore(std::istream& fin){
		std::string s; fin >> s; // discard version number 
	}

};


#endif
