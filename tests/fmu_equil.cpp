#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

#include <solver.h>
#include "pspm_interface.h"

int main(){

	PSPM_Plant p1;
	p1.initParamsFromFile("tests/params/p.ini");
	p1.set_size(0.01);

	PSPM_Dynamic_Environment E;
	E.metFile = "tests/data/MetData_AmzFACE_Monthly_2000_2015_PlantFATE.csv";
	E.co2File = "tests/data/CO2_AMB_AmzFACE2000_2100.csv";
	E.init();
	E.print(0);

	Species<PSPM_Plant> spp(p1);

	Solver S(SOLVER_FMU);
	S.setEnvironment(&E);
	S.addSpecies(30, 0.01, 1.5, true, &spp, 3, 1);
	//S.get_species(0)->set_bfin_is_u0in(true);	// say that input_birth_flux is u0
	S.resetState();
	S.initialize();
	S.print();
	
//	ofstream fout("fmu_Redmodel.txt");

//	for (double t=0.05; t <= 5000; t=t+100) {
//		S.step_to(t);
//		fout << S.current_time << "\t" << S.newborns_out()[0] << "\t";
//		//cout << S.current_time << " " [><< S.u0_out()<] << "\n";
//		for (auto y : S.state) fout << y << "\t";
//		fout << endl;
//	}

//	fout.close();

//	// Expected 44.3530812 
//	cout << S.newborns_out()[0] << endl; 
//	//if (abs(S.u0_out()[0] - 1.468232) < 1e-5) return 0;
//	//else return 1;

}

