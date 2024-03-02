#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include "solver.h"
using namespace std;

#include "red.h"

int main(){

	Species<RED_Plant> spp;
	LightEnvironment E;

	Solver S(SOLVER_IFMU);
	S.control.ifmu_order = 1;
	S.control.ode_ifmu_stepsize = 1;

	S.setEnvironment(&E);
	S.addSpecies({150}, {1}, {1e6}, {true}, &spp, 0);
	//S.get_species(0)->set_bfin_is_u0in(true);	// say that input_birth_flux is u0
	S.initialize();
	//S.print();
	
	ofstream fout("ifmu_Redmodel.txt");

	for (double t=0.05; t <= 5000; t=t+10) {
		S.step_to(t);
		fout << S.current_time << "\t" << S.newborns_out(t)[0] << "\t";
		cout << S.current_time << "\t" << S.newborns_out(t)[0] << "\n";
		
		for (auto y : S.state) fout << y << "\t";
		fout << endl;
	}

	fout.close();

	// Expected 44.3530812 
	cout << S.newborns_out(5000)[0] << endl; 
	//if (abs(S.u0_out()[0] - 1.468232) < 1e-5) return 0;
	//else return 1;

}

