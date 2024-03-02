#ifndef  PSPM_PSPM_ENVIRONMENT_H_
#define  PSPM_PSPM_ENVIRONMENT_H_

class Solver;

class EnvironmentBase{
	public:
	virtual ~EnvironmentBase(){};
	virtual void computeEnv(double t, Solver * sol, std::vector<double>::iterator S, std::vector<double>::iterator dSdt) = 0;
/*	virtual void calcRatesSystem(double t, vector<double>::iterator S, vector<double>::iterator dSdt) = 0;*/
};


#endif

