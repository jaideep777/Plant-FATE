#include <iostream>
#include <vector>
#include <cassert>
#include <environment_base.h>
#include <individual_base.h>
using namespace std;


class LightEnv : public EnvironmentBase{
	public:
	double E = 0.95;

	void computeEnv(double t, Solver * S, vector<double>::iterator s, vector<double>::iterator dsdt){
		E = 0.95 + 0.05*t;
	}

};



class Plant : public IndividualBase<1>{
	public:
	double lma = 30;

	double height = -99;
	double crown_area = -99;
	double root_mass = -99;

	vector<string> varnames = {"lma|", "ht", "cr", "root"};

	Plant(){
		lma = 10;
	}

	Plant(double _lma){
		lma = _lma;
	}
	
	void set_size(const array<double,1>& x){
		height = x[0];
		crown_area = height*height;
	}

	void preCompute(double t, void * env){
	}

	array<double,1> growthRate(double t, void * env){
		cout << "Plant::g(): " << x[0] << " " << t << " " << ((LightEnv*)env)->E << "\n";
		return {x[0]*((LightEnv*)env)->E};
	}
	double mortalityRate(double t, void * env){
		return -0.5;
	}
	double birthRate(double t, void * env){
		return 1;
	}
	
	double establishmentProbability(double t, void  * _env){
		return 1;
	}

	double init_density(void * _env, double bf){
		return 5/(x[0]+0.5);
	}

	void init_accumulators(double t, void * env){
		root_mass = 0.03 + 2*((LightEnv*)env)->E;
	}

	vector<double>::iterator set_accumulators(vector<double>::iterator &it){
		root_mass = *it++;
		return it;
	}

	vector<double>::iterator get_accumulators(vector<double>::iterator &it){
		*it++ = root_mass;
		return it;
	}

	vector<double>::iterator get_accumulatorRates(vector<double>::iterator &it){
		*it++ = -0.1*(root_mass-1);
		return it;
	}

	void print(std::ostream& out = std::cout) const {
		out << std::setw(6) << std::setprecision(4) 
		    << lma << "\t" << height << "\t" << crown_area << "\t" << root_mass << "\t";
	}
	
	void save(std::ostream& fout){
		fout << "TestPlant::v1" << "   ";
		fout << lma << "   ";
	}
	void restore(std::istream& fin){
		std::string s; fin >> s; // discard version number
		assert(s == "TestPlant::v1");
		fin >> lma;
	}
	
};


class Insect : public IndividualBase<2>{
	public:
	double w; // weight
	double e; // energy
	double wingspan = 10;
	
	double m = -1, f = -1;
	array<double,2> g;

	vector<string> varnames = {"wing", "f"};
	
	double init_density(void * _env, double bf){
		return x[0]*x[1]/10;
	}

	void set_size(const array<double,2>& x){
		w = x[0];
		e = x[1];
		wingspan = 2*w;	
	}

	void preCompute(double t, void * env){
		cout << "Insect::Precompute()\n";
		LightEnv * le = static_cast<LightEnv*>(env);
		f = 0.1*w;
		g = {w*e*0.1, le->E*w*0.1};
		m = w/e;
	}

	array<double,2> growthRate(double t, void * env){
		cout << "Insect::g()\n";
		return g;
	}
	double mortalityRate(double t, void * env){
		cout << "Insect::m()\n";
		return m;
	}
	double birthRate(double t, void * env){
		cout << "Insect::f()\n";
		return f;
	}

	double establishmentProbability(double t, void  * _env){
		return 1;
	}

	void print(std::ostream& out = std::cout) const {
		out << std::setw(12) << std::setprecision(4) 
		    << w << '\t' << e << '\t'
			<< wingspan << '\t' << g << '\t' << m << '\t' << f << "\t";
	}
};
