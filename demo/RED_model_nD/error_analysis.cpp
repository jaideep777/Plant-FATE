#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono>
#include "solver.h"
using namespace std;

#include "red.h"

inline std::vector <double> seq(double from, double to, int len){
	std::vector<double> x(len);
	for (size_t i=0; i<len; ++i) x[i] = from + i*(to-from)/(len-1);
	return x;
}


int main(){

    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;

	{
	cout << "running EBT...\n";
	// EBT
	ofstream ferr("ebt_error_analysis.txt");
	ferr << "N0\tNf\tdt\tB\tEb\ttsys\n";
	
	for (int i=3; i<10; ++i){
		int N0 = pow(2,i);
		double Dt = 10*pow(2, 5-i);
		cout << "N0 = " << N0 << ", Dt = " << Dt << endl;
	
		Species<RED_Plant> spp;
		LightEnvironment E;

		Solver S(SOLVER_EBT);
		S.control.ebt_ucut = 1e-20;

		S.setEnvironment(&E);
		S.addSpecies({100}, {1}, {1e6}, {true}, &spp, 0);
		//S.get_species(0)->set_bfin_is_u0in(true);	// say that input_birth_flux is u0
		S.resetState();
		S.initialize();
		//S.print();

		
		auto t1 = high_resolution_clock::now();
		for (double t=0; t <= 5000; t=t+Dt) {
			S.step_to(t);
		}
	    auto t2 = high_resolution_clock::now();
   		duration<double, std::milli> ms_double = t2 - t1;
   
		double B = S.integrate_x([&S](int i, double t){return S.species_vec[0]->getX(i);}, S.current_time, 0);
		ferr << N0 << "\t" << S.species_vec[0]->xsize() << "\t" << Dt << "\t" << B << "\t" << fabs(B-13933.16)/13933.16 << "\t" << ms_double.count() << "\n";
	}

	ferr.close();
	}
	

	{
	cout << "running IEBT...\n";
	// IEBT
	ofstream ferr("iebt_error_analysis.txt");
	ferr << "N0\tNf\tdt\tB\tEb\ttsys\n";
	
	for (int i=3; i<10; ++i){
		int N0 = pow(2,i);
		double Dt = 10*pow(2, 5-i);
		cout << "N0 = " << N0 << ", Dt = " << Dt << endl;
	
		Species<RED_Plant> spp;
		LightEnvironment E;

		Solver S(SOLVER_IEBT);
		S.control.ebt_ucut = 1e-20;

		S.setEnvironment(&E);
		S.addSpecies({100}, {1}, {1e6}, {true}, &spp, 0);
		//S.get_species(0)->set_bfin_is_u0in(true);	// say that input_birth_flux is u0
		S.resetState();
		S.initialize();
		//S.print();

		
		auto t1 = high_resolution_clock::now();
		for (double t=0; t <= 5000; t=t+Dt) {
			S.step_to(t);
		}
	    auto t2 = high_resolution_clock::now();
   		duration<double, std::milli> ms_double = t2 - t1;
   
		double B = S.integrate_x([&S](int i, double t){return S.species_vec[0]->getX(i);}, S.current_time, 0);
		ferr << N0 << "\t" << S.species_vec[0]->xsize() << "\t" << Dt << "\t" << B << "\t" << fabs(B-13933.16)/13933.16 << "\t" << ms_double.count() << "\n";
	}

	ferr.close();
	}

	
	
	{
	cout << "running FMU...\n";
	// FMU
	ofstream ferr("fmu_error_analysis.txt");
	ferr << "N0\tNf\tdt\tB\tEb\ttsys\n";
	
	for (int i=3; i<11; ++i){
		int N0 = pow(2,i);
		cout << "N0 = " << N0 << endl;
			
		Species<RED_Plant> spp;
		LightEnvironment E;

		Solver S(SOLVER_FMU);
		S.setEnvironment(&E);
		S.addSpecies(N0, 1, 1e6, true, &spp, 0);
		//S.get_species(0)->set_bfin_is_u0in(true);	// say that input_birth_flux is u0
		S.resetState();
		S.initialize();
		//S.print();
			
		auto t1 = high_resolution_clock::now();
		for (double t=0; t <= 5000; t=t+10) {
			S.step_to(t);
		}
	    auto t2 = high_resolution_clock::now();
   		duration<double, std::milli> ms_double = t2 - t1;
		
		double B = S.integrate_x([&S](int i, double t){return S.species_vec[0]->getX(i);}, S.current_time, 0);
		ferr << N0 << "\t" << S.species_vec[0]->xsize() << "\t" << 0 << "\t" << B << "\t" << fabs(B-13933.16)/13933.16 << "\t" << ms_double.count() << "\n";
	}

	ferr.close();
	}



	{
	cout << "running IFMU...\n";
	// IFMU
	ofstream ferr("ifmu_error_analysis.txt");
	ferr << "N0\tNf\tdt\tB\tEb\ttsys\n";
	
		for (int i=3; i<11; ++i){
		int N0 = pow(2,i);
		cout << "N0 = " << N0 << endl;
		
		Species<RED_Plant> spp;
		LightEnvironment E;

		Solver S(SOLVER_IFMU);
		S.control.ifmu_order = 1;
		S.control.ode_ifmu_stepsize = 0.1;

		S.setEnvironment(&E);
		S.addSpecies(N0, 1, 1e6, true, &spp, 0);
		//S.get_species(0)->set_bfin_is_u0in(true);	// say that input_birth_flux is u0
		S.resetState();
		S.initialize();
		//S.print();
			
		auto t1 = high_resolution_clock::now();
		for (double t=0; t <= 5000; t=t+10) {
			S.step_to(t);
		}
	    auto t2 = high_resolution_clock::now();
   		duration<double, std::milli> ms_double = t2 - t1;
		
		
		double B = S.integrate_x([&S](int i, double t){return S.species_vec[0]->getX(i);}, S.current_time, 0);
		ferr << N0 << "\t" << S.species_vec[0]->xsize() << "\t" << 0 << "\t" << B << "\t" << fabs(B-13933.16)/13933.16 << "\t" << ms_double.count() << "\n";
	}

	ferr.close();
	}


	{
	cout << "running IFMU(O2)...\n";
	// IFMU
	ofstream ferr("ifmu2_error_analysis.txt");
	ferr << "N0\tNf\tdt\tB\tEb\ttsys\n";
	
		for (int i=3; i<11; ++i){
		int N0 = pow(2,i);
		cout << "N0 = " << N0 << endl;
		
		Species<RED_Plant> spp;
		LightEnvironment E;

		Solver S(SOLVER_IFMU);
		S.control.ifmu_order = 2;
		S.control.ode_ifmu_stepsize = 0.1;

		S.setEnvironment(&E);
		S.addSpecies(N0, 1, 1e6, true, &spp, 0);
		//S.get_species(0)->set_bfin_is_u0in(true);	// say that input_birth_flux is u0
		S.resetState();
		S.initialize();
		//S.print();
			
		auto t1 = high_resolution_clock::now();
		for (double t=0; t <= 5000; t=t+10) {
			S.step_to(t);
		}
	    auto t2 = high_resolution_clock::now();
   		duration<double, std::milli> ms_double = t2 - t1;
		
		
		double B = S.integrate_x([&S](int i, double t){return S.species_vec[0]->getX(i);}, S.current_time, 0);
		ferr << N0 << "\t" << S.species_vec[0]->xsize() << "\t" << 0 << "\t" << B << "\t" << fabs(B-13933.16)/13933.16 << "\t" << ms_double.count() << "\n";
	}

	ferr.close();
	}


	{
	cout << "running CM...\n";
	// EBT
	ofstream ferr("cm_error_analysis.txt");
	ferr << "N0\tNf\tdt\tB\tEb\ttsys\n";
	
	for (int i=3; i<10; ++i){
		int N0 = pow(2,i);
		double Dt = 10*pow(2, 5-i);
		cout << "N0 = " << N0 << ", Dt = " << Dt << endl;
	
		Species<RED_Plant> spp;
		LightEnvironment E;

		Solver S(SOLVER_CM);
		S.setEnvironment(&E);
		S.addSpecies(N0, 1, 1e6, true, &spp, 0);
		//S.get_species(0)->set_bfin_is_u0in(true);	// say that input_birth_flux is u0
		//S.control.max_cohorts = N0;
		S.resetState();
		S.initialize();
		//S.print();

		
		auto t1 = high_resolution_clock::now();
		for (double t=0; t <= 5000; t=t+Dt) {
			S.step_to(t);
			// cout << S.current_time << "\t" << S.newborns_out(t)[0] << "\t";
			// cout << S.current_time << " " << S.species_vec[0]->xsize() << "\n";
		}
	    auto t2 = high_resolution_clock::now();
    	duration<double, std::milli> ms_double = t2 - t1;
    
		double B = S.integrate_x([&S](int i, double t){return S.species_vec[0]->getX(i);}, S.current_time, 0);
		ferr << N0 << "\t" << S.species_vec[0]->xsize() << "\t" << Dt << "\t" << B << "\t" << fabs(B-13933.16)/13933.16 << "\t" << ms_double.count() << "\n";
	}

	ferr.close();
	}

	{
	cout << "running ICM...\n";
	// EBT
	ofstream ferr("icm_error_analysis.txt");
	ferr << "N0\tNf\tdt\tB\tEb\ttsys\n";
	
	for (int i=3; i<10; ++i){
		int N0 = pow(2,i);
		double Dt = 10*pow(2, 5-i);
		cout << "N0 = " << N0 << ", Dt = " << Dt << endl;
	
		Species<RED_Plant> spp;
		LightEnvironment E;

		Solver S(SOLVER_ICM);
		S.control.cm_use_log_densities = false;

		S.setEnvironment(&E);
		S.addSpecies(N0, 1, 1e6, true, &spp, 0);
		//S.get_species(0)->set_bfin_is_u0in(true);	// say that input_birth_flux is u0
		//S.control.max_cohorts = N0;
		S.resetState();
		S.initialize();
		//S.print();

		
		auto t1 = high_resolution_clock::now();
		for (double t=0; t <= 5000; t=t+Dt) {
			S.step_to(t);
			// cout << S.current_time << "\t" << S.newborns_out(t)[0] << "\t";
			// cout << S.current_time << " " << S.species_vec[0]->xsize() << "\n";
		}
	    auto t2 = high_resolution_clock::now();
    	duration<double, std::milli> ms_double = t2 - t1;
    
		double B = S.integrate_x([&S](int i, double t){return S.species_vec[0]->getX(i);}, S.current_time, 0);
		ferr << N0 << "\t" << S.species_vec[0]->xsize() << "\t" << Dt << "\t" << B << "\t" << fabs(B-13933.16)/13933.16 << "\t" << ms_double.count() << "\n";
	}

	ferr.close();
	}

	{
	cout << "running ABM...\n";
	// ABM
	ofstream ferr("abm_error_analysis.txt");
	ferr << "N0\tNf\tdt\tB\tEb\ttsys\n";
	
	for (int i=3; i<11; ++i){
		int N0 = pow(2,i);
		cout << "N0 = " << N0 << endl;
		Species<RED_Plant> spp;
		LightEnvironment E;

		Solver S(SOLVER_ABM);
		S.control.abm_n0 = N0;
		S.control.abm_stepsize = 0.5;

		S.setEnvironment(&E);
		S.addSpecies({100}, {1}, {1e6}, {true}, &spp, 0);
		//S.get_species(0)->set_bfin_is_u0in(true);	// say that input_birth_flux is u0
		S.resetState();
		S.initialize();
		//S.print();
			
		auto t1 = high_resolution_clock::now();
		for (double t=0; t <= 5000; t=t+10) {
			S.step_to(t);
		}
		auto t2 = high_resolution_clock::now();
   		duration<double, std::milli> ms_double = t2 - t1;
	
		
		double B = S.integrate_x([&S](int i, double t){return S.species_vec[0]->getX(i);}, S.current_time, 0);
		ferr << N0 << "\t" << S.species_vec[0]->xsize() << "\t" << 0 << "\t" << B << "\t" << fabs(B-13933.16)/13933.16 << "\t" << ms_double.count() << "\n";
	}

	ferr.close();
	}	

}


