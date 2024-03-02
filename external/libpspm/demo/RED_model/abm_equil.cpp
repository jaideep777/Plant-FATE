#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include "solver.h"
using namespace std;

#include "red.h"

inline std::vector <double> logseq(double from, double to, int len){
	std::vector<double> x(len);
	for (size_t i=0; i<len; ++i) x[i] = exp(log(from) + i*(log(to)-log(from))/(len-1));
	return x;
}


int main(){

	Species<RED_Plant> spp;
	LightEnvironment E;

	Solver S(SOLVER_ABM);
	S.control.abm_stepsize = 1;
	S.control.abm_n0 = 1000;
	S.setEnvironment(&E);

	S.addSpecies({100}, {1}, {1e6}, {true}, &spp, 0);
	//S.get_species(0)->set_bfin_is_u0in(true);	// say that input_birth_flux is u0
	S.initialize();
	//S.print();
	
	
	ofstream fout("abm_Redmodel.txt");

	for (double t=0; t <= 5000; t=t+10) {
		S.step_to(t);
		fout << S.current_time << "\t" << S.newborns_out(t)[0] << "\t";
		cout << S.current_time << " " << S.species_vec[0]->xsize() << "\n";
		//cout << S.current_time << " " [><< S.u0_out()<] << "\n";
		
		vector <double> dist = S.getDensitySpecies1D(0, 0, logseq(1, 1e6, 150));
		for (auto y : dist) fout << y << "\t";
		fout << endl;
	}

	fout.close();

	// expected 38.1128953 (numerical R), 37.5845 (analytical)
	cout << S.newborns_out(5000)[0] << endl; 
	//if (abs(S.u0_out()[0] - 1.468232) < 1e-5) return 0;
	//else return 1;

}

