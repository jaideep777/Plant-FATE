#include <iomanip>
#include <fstream>
#include <random>

#include "life_history.h"
using namespace std;

int main(int argc, char ** argv){

	if (argc < 2){
		cout << "Plant-FATE error: pf needs at least 1 argument. Syntax: ./pfate <params_file.ini> [<n_years>]\n";
		return 1;
	}

	string pfile = argv[1];
	
	double nyears = 500;
	if (argc > 2) nyears = stod(argv[2]);

	if (nyears <=0 ){
		cout << "Plant-FATE error: number of years should be > 0\n";
	}

	cout << setprecision(12);
	
	pfate::LifeHistoryOptimizer lho(pfile);
	// lho.C.init_co2(414);
	lho.init();
	lho.C.Climate::print(0);
	double total_prod = lho.P.get_biomass();
	cout << "Starting biomass = " << total_prod << "\n";
	cout << "Mortality until seedling stage = " << lho.P.state.mortality << "\n";

	ofstream fout("assim1.txt");
	double dt = 0.1;
	lho.printHeader(fout);
	for (double t=2000; t<=2000+nyears; t=t+dt){
		lho.grow_for_dt(t, dt);
		lho.printState(t+dt, fout);
		// lho.C.Climate::print(t);
	}
	fout.close();
	// fswp.close();

	cout << "At last t: " << "\n" 
		 << "  Total biomass    = " << lho.P.get_biomass() << "\n"
		 << "  Total litter     = " << lho.litter_pool << "\n"
		 << "  Total reproduc   = " << lho.rep << "\n"
		 << "  Total bio+lit+rep = " << lho.P.get_biomass() + lho.litter_pool + lho.rep << "\n"
		 << "  Total production = " << lho.prod << "\n";
	

	double rel_error = abs((lho.P.get_biomass()+lho.litter_pool+lho.rep)/lho.prod - 1);
	cout << "Relative error in biomass accounting = " << rel_error << endl;
	
	double fitness1 = lho.seeds;
	cout << "Fitness = " << fitness1 << '\n';

	return 0;
}
