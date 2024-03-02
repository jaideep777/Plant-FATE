#include <iostream>
#include <cmath>
#include <environment_base.h>
#include <individual_base.h>
using namespace std;



class WaveEnv : public EnvironmentBase{
	public:
	void computeEnv(double t, Solver * S, vector<double>::iterator s, vector<double>::iterator dsdt){
	}

};

class Wave : public IndividualBase<2>{
	public:
	double init_density(void * _env, double bf){
		return exp(-10*(x[0]-2)*(x[0]-2) - 10*(x[1]-4)*(x[1]-4));
		// return (x[0] > 2 && x[0] < 3 && x[1] > 3 && x[1] < 4)? 1 : 0;
	}

	void set_size(const array<double,2>& x){
	}

	array<double,2> growthRate(double t, void * env){
		return {1,2};
	}
	double mortalityRate(double t, void * env){
		return 0;
	}
	double birthRate(double t, void * env){
		return 0.5;
	}

};

