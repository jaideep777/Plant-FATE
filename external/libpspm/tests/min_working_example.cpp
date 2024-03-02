#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>

#include <individual_base.h>
#include <solver.h>

using namespace std;

class Environment : public EnvironmentBase{
	public:
	double E = 0;

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



class Individual : public IndividualBase<1>{
	public:
	double init_density(void * _env, double bf){
		double _x = x[0];
		return pow(1-_x,2)/pow(1+_x,4) + (1-_x)/pow(1+_x,3);
	}

	std::array<double,1> growthRate(double t, void * _env){
		Environment* env = (Environment*)_env;
		double _x = x[0];
		double E = env->E;
		double a = 0.16+0.22*exp(-0.225*t*t);
		return {0.225*(1-_x*_x)*(E/(1+E*E))*t*(1+a*a)/a};
	}

	double mortalityRate(double t, void * _env){
		Environment* env = (Environment*)_env;
		double _x = x[0];
		double E = env->E;
		double a = 0.16+0.22*exp(-0.225*t*t);
		return 1.35*t*E/a;
	}

	double birthRate(double t, void * _env){
		Environment* env = (Environment*)_env;
		double _x = x[0];
		double E = env->E;
		double oneplusa = 1.16+0.22*exp(-0.225*t*t);
		double a = 0.16+0.22*exp(-0.225*t*t);
		double n1 = 0.225*t*_x*_x*(1-_x)*(1-_x)*E/(1+E)/(1+E)*oneplusa*oneplusa/a;
		double n2 = (1+exp(-0.225*t*t))/(61-88*log(2)+(38*log(2)-79.0/3)*exp(-0.225*t*t));
		return n1*n2;
	}


	void set_size(const std::array <double, 1>& _x){
	}


};


int main(){
	Environment E;
	Species<Individual> spp;

	Solver S(SOLVER_FMU, "rk45ck");
	S.setEnvironment(&E);
	S.addSpecies({25}, {0}, {1}, {false}, &spp, 0, -1);

	for (double t=0.05; t <= 8; t=t+0.05) {
		S.step_to(t);
		cout << S.current_time << " " << S.u0_out(t)[0] << "\n";
	}

	if (abs(S.u0_out(S.current_time)[0] - 0.958400) < 1e-6) return 0;
	else return 1;

}


