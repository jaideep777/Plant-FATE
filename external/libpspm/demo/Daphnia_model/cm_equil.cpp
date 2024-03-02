#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include "solver.h"
using namespace std;

#include "daphnia.h"

inline std::vector <double> seq(double from, double to, int len){
	std::vector<double> x(len);
	for (size_t i=0; i<len; ++i) x[i] = from + i*(to-from)/(len-1);
	return x;
}


int main(){

	Species<Daphnia> spp;
	Environment E;

	Solver S(SOLVER_CM);
	S.control.max_cohorts = 200;
	S.control.update_cohorts = true;
	S.control.cm_remove_cohorts = true;
	S.control.cohort_insertion_dt = 0.5;
	// S.control.sync_cohort_insertion = false;
	S.control.cm_use_log_densities = true;

	S.setEnvironment(&E);

	S.addSystemVariables({E.K});  // this can be done either before or after addSpecies()
	S.addSpecies({100}, {0}, {1}, {false}, &spp, 0, -1);

	S.initialize();
	//S.print();
	
	
	ofstream fout("cm_Daphnia.txt");

	for (double t=0.05; t <= 100; t=t+0.5) {
		S.step_to(t);
		fout << S.current_time << "\t" << S.newborns_out(t)[0] << "\t" << E.S << "\t";
		cout << S.current_time << " " << S.state[0] << " " << S.species_vec[0]->xsize() << "\n";
		
		vector <double> dist = S.getDensitySpecies1D(0, 0, seq(0,1,300));
		for (auto y : dist) fout << y << "\t";
		
		fout << endl;
	}

	fout.close();

	// expected 38.1128953 (numerical R), 37.5845 (analytical)
	cout << S.newborns_out(100)[0] << endl; 
	//if (abs(S.u0_out()[0] - 1.468232) < 1e-5) return 0;
	//else return 1;
	// S.print();
}

