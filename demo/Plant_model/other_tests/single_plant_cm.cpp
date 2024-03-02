#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
using namespace std;

#include <solver.h>
#include "pspm_environment.h"
#include "pspm_plant.h"


vector<double> generateDefaultCohortSchedule(double max_time){

	vector<double> tvec;

	const double multiplier=0.2, min_step_size=1e-5, max_step_size=2;
	
	assert(min_step_size > 0 && "The minimum step size must be greater than zero");
	
	double dt = 0.0, time = 0.0;
	tvec.push_back(time);
	while (time <= max_time) {
		dt = exp2(floor(log2(time * multiplier)));
		time += min(max(dt, min_step_size), max_step_size);
		tvec.push_back(time);
	}

	// Drop the last time; that's not going to be needed:
	if (tvec.size() >=1) 	// JAI: added to avoid overflow warning
		tvec.resize(tvec.size() - 1);

	return tvec;
}



int main(){
	
	//initPlantParameters(plant::par);
	
	FixedEnvironment env(1);	
	
	PSPM_Plant p;
	p.lma = 0.1978791;
	p.initParameters();
	p.vars.height = p.par.height_0; //0.3257146; //0.3920458; //0.3441948;
	p.vars.area_leaf = p.par.area_leaf_0; 

	cout << p << endl;

    Species<PSPM_Plant> s1(p);
    //M.p = M.seed = p;
	s1.print(); 
	
	Solver S(SOLVER_CM);
	S.control.cm_use_log_densities = true;
	S.control.ode_eps = 1e-4;
	S.control.update_cohorts = false;
	S.setEnvironment(&env);

	S.addSpecies(vector<double>(1, p.vars.height), &s1, 4, 1);
	
	S.resetState();
	S.initialize();

	S.print();

	vector <double> times = generateDefaultCohortSchedule(105.32);
	for (auto t : times) cout << t << " "; cout << endl;

	ofstream fout("out_single_plant_cm.txt");
	
	for (size_t i=0; i < times.size(); ++i){

		S.step_to(times[i]);		
		fout << times[i] << "\t"; ((Species<PSPM_Plant>*)S.species_vec[0])->cohorts[0].print(fout); fout << "\n";

	}

	S.print();

	fout.close();

}

