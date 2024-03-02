#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono>
#include "solver.h"
using namespace std;

#include "daphnia.h"

inline std::vector <double> seq(double from, double to, int len){
	std::vector<double> x(len);
	for (size_t i=0; i<len; ++i) x[i] = from + i*(to-from)/(len-1);
	return x;
}

    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;

void run_ebt(){
	cout << "running EBT...\n";
	// EBT
	ofstream ferr("ebt_error_analysis.txt");
	ferr << "N0\tNf\tdt\tB\tEb\ttsys\n";
	
	for (int i=3; i<12; ++i){
		int N0 = pow(2,i);
		double Dt = 0.2*pow(2, 9-i);
		cout << "N0 = " << N0 << ", Dt = " << Dt << endl;

		Species<Daphnia> spp;
		Environment E;

		Solver S(SOLVER_EBT);
		S.control.cohort_insertion_dt = Dt;
		S.control.ebt_ucut = 1e-10;

		S.setEnvironment(&E);
		S.addSpecies({100}, {0}, {1}, {false}, &spp, 0, -1);
		S.addSystemVariables({E.K});  // this can be done either before or after addSpecies()

		S.initialize();
		//S.print();
		
		auto t1 = high_resolution_clock::now();
		S.step_to(150);
		auto t2 = high_resolution_clock::now();
		duration<double, std::milli> ms_double = t2 - t1;
   
		double B = S.state_integral([&S](int i, double t){return S.species_vec[0]->getX(i)[0];}, S.current_time, 0);
		ferr << N0 << "\t" << S.species_vec[0]->xsize() << "\t" << Dt << "\t" << B << "\t" << fabs(B-1.298077)/1.298077 << "\t" << ms_double.count() << "\n";
	}

	ferr.close();
}
	
void run_iebt(){
	cout << "running IEBT...\n";
	// EBT
	ofstream ferr("iebt_error_analysis.txt");
	ferr << "N0\tNf\tdt\tB\tEb\ttsys\n";
	
	for (int i=3; i<12; ++i){
		int N0 = pow(2,i);
		double Dt = 0.2*pow(2, 9-i);
		cout << "N0 = " << N0 << ", Dt = " << Dt << endl;

		Species<Daphnia> spp;
		Environment E;

		Solver S(SOLVER_IEBT);
		S.control.cohort_insertion_dt = Dt;
		S.control.ebt_ucut = 1e-10;

		S.setEnvironment(&E);
		S.addSpecies({100}, {0}, {1}, {false}, &spp, 0, -1);
		S.addSystemVariables({E.K});  // this can be done either before or after addSpecies()

		S.initialize();
		//S.print();
		
		auto t1 = high_resolution_clock::now();
		S.step_to(150);
		auto t2 = high_resolution_clock::now();
		duration<double, std::milli> ms_double = t2 - t1;
   
		double B = S.state_integral([&S](int i, double t){return S.species_vec[0]->getX(i)[0];}, S.current_time, 0);
		ferr << N0 << "\t" << S.species_vec[0]->xsize() << "\t" << Dt << "\t" << B << "\t" << fabs(B-1.298077)/1.298077 << "\t" << ms_double.count() << "\n";
	}

	ferr.close();
}
		
	
void run_fmu(){
	cout << "running FMU...\n";
	// FMU
	ofstream ferr("fmu_error_analysis.txt");
	ferr << "N0\tNf\tdt\tB\tEb\ttsys\n";
	
	for (int i=3; i<11; ++i){
		int N0 = pow(2,i);
		cout << "N0 = " << N0 << endl;
	
		Species<Daphnia> spp;
		Environment E;

		Solver S(SOLVER_FMU);
		S.setEnvironment(&E);
		S.addSpecies({N0}, {0}, {1}, {false}, &spp, 0, -1);
		S.addSystemVariables({E.K});  // this can be done either before or after addSpecies()

		S.initialize();
		//S.print();
		
		auto t1 = high_resolution_clock::now();
		S.step_to(150);
	    auto t2 = high_resolution_clock::now();
		duration<double, std::milli> ms_double = t2 - t1;
		
		double B = S.state_integral([&S](int i, double t){return S.species_vec[0]->getX(i)[0];}, S.current_time, 0);
		ferr << N0 << "\t" << S.species_vec[0]->xsize() << "\t" << 0 << "\t" << B << "\t" << fabs(B-1.298077)/1.298077 << "\t" << ms_double.count() << "\n";
	}

	ferr.close();
}



void run_ifmu(){
	cout << "running IFMU...\n";
	// IFMU
	ofstream ferr("ifmu_error_analysis.txt");
	ferr << "N0\tNf\tdt\tB\tEb\ttsys\n";
	
		for (int i=3; i<13; ++i){
		int N0 = pow(2,i);
		cout << "N0 = " << N0 << endl;
		
		Species<Daphnia> spp;
		Environment E;

		Solver S(SOLVER_IFMU);
		S.control.ode_ifmu_stepsize = 0.02;
		S.control.ifmu_order = 1;
		
		S.setEnvironment(&E);
		S.addSpecies({N0}, {0}, {1}, {false}, &spp, 0, -1);
		S.addSystemVariables({E.K});  // this can be done either before or after addSpecies()

		S.initialize();
		//S.print();
		
		auto t1 = high_resolution_clock::now();
		S.step_to(150);
	    auto t2 = high_resolution_clock::now();
    	duration<double, std::milli> ms_double = t2 - t1;
		
		
		double B = S.state_integral([&S](int i, double t){return S.species_vec[0]->getX(i)[0];}, S.current_time, 0);
		ferr << N0 << "\t" << S.species_vec[0]->xsize() << "\t" << 0 << "\t" << B << "\t" << fabs(B-1.298077)/1.298077 << "\t" << ms_double.count() << "\n";
	}

	ferr.close();
}


	// {
	// cout << "running IFMU(O2)...\n";
	// // IFMU
	// ofstream ferr("ifmu2_error_analysis.txt");
	// ferr << "N0\tNf\tdt\tB\tEb\ttsys\n";
	
	// for (int i=3; i<13; ++i){
	// 	int N0 = pow(2,i);
	// 	cout << "N0 = " << N0 << endl;
		
	// 	Species<Daphnia> spp;
	// 	Environment E;

	// 	Solver S(SOLVER_IFMU);
	// 	S.control.ode_ifmu_stepsize = 0.02;
	// 	S.control.ifmu_order = 2;
		
	// 	S.setEnvironment(&E);
	// 	S.addSpecies({N0}, {0}, {1}, {false}, &spp, 0, -1);
	// 	S.addSystemVariables({E.K});  // this can be done either before or after addSpecies()

	// 	S.initialize();
	// 	//S.print();
		
	// 	auto t1 = high_resolution_clock::now();
	// 	S.step_to(150);
	//     auto t2 = high_resolution_clock::now();
    // 	duration<double, std::milli> ms_double = t2 - t1;
		
		
	// 	double B = S.state_integral([&S](int i, double t){return S.species_vec[0]->getX(i)[0];}, S.current_time, 0);
	// 	ferr << N0 << "\t" << S.species_vec[0]->xsize() << "\t" << 0 << "\t" << B << "\t" << fabs(B-1.298077)/1.298077 << "\t" << ms_double.count() << "\n";
	// }

	// ferr.close();
	// }


void run_cm(){
	cout << "running CM...\n";
	// IFMU
	ofstream ferr("cm_error_analysis.txt");
	ferr << "N0\tNf\tdt\tB\tEb\ttsys\n";
	
	for (int i=3; i<12; ++i){
		int N0 = pow(2,i);
		double Dt = 0.2*pow(2, 9-i);
		cout << "N0 = " << N0 << ", Dt = " << Dt << endl;
		
		Species<Daphnia> spp;
		Environment E;

		Solver S(SOLVER_CM);
		S.control.cohort_insertion_dt = Dt;
		S.control.max_cohorts = 1000;
		S.control.cm_remove_cohorts = true;
		S.control.ebt_ucut = 1e-10;
		S.control.cm_dxcut = 1e-5;

		S.setEnvironment(&E);
		S.addSpecies({100}, {0}, {1}, {false}, &spp, 0, -1);
		S.addSystemVariables({E.K});  // this can be done either before or after addSpecies()

		S.initialize();
		//S.print();
	
		auto t1 = high_resolution_clock::now();
		for (double t=0.05; t <= 150; t=t+Dt) {
			S.step_to(t);
		}
	    auto t2 = high_resolution_clock::now();
    	duration<double, std::milli> ms_double = t2 - t1;
		
		
		double B = S.state_integral([&S](int i, double t){return S.species_vec[0]->getX(i)[0];}, S.current_time, 0);
		ferr << N0 << "\t" << S.species_vec[0]->xsize() << "\t" << 0 << "\t" << B << "\t" << fabs(B-1.298077)/1.298077 << "\t" << ms_double.count() << "\n";
	}

	ferr.close();
}

void run_icm(){
	cout << "running ICM...\n";
	// IFMU
	ofstream ferr("icm_error_analysis.txt");
	ferr << "N0\tNf\tdt\tB\tEb\ttsys\n";
	
	for (int i=3; i<12; ++i){
		int N0 = pow(2,i);
		double Dt = 0.2*pow(2, 9-i);
		cout << "N0 = " << N0 << ", Dt = " << Dt << endl;
		
		Species<Daphnia> spp;
		Environment E;

		Solver S(SOLVER_ICM);
		S.control.cm_use_log_densities = false;
		S.control.max_cohorts = 1000;
		S.control.cm_remove_cohorts = true;
		S.control.ebt_ucut = 1e-10;
		S.control.cm_dxcut = 1e-5;
		S.control.cohort_insertion_dt = Dt;

		S.setEnvironment(&E);
		S.addSpecies({100}, {0}, {1}, {false}, &spp, 0, -1);
		S.addSystemVariables({E.K});  // this can be done either before or after addSpecies()

		S.initialize();
		//S.print();
	
		auto t1 = high_resolution_clock::now();
		for (double t=0.05; t <= 150; t=t+Dt) {
			S.step_to(t);
		}
	    auto t2 = high_resolution_clock::now();
    	duration<double, std::milli> ms_double = t2 - t1;
		
		
		double B = S.state_integral([&S](int i, double t){return S.species_vec[0]->getX(i)[0];}, S.current_time, 0);
		ferr << N0 << "\t" << S.species_vec[0]->xsize() << "\t" << 0 << "\t" << B << "\t" << fabs(B-1.298077)/1.298077 << "\t" << ms_double.count() << "\n";
	}

	ferr.close();
}


void run_abm(){
	cout << "running ABM...\n";
	// ABM
	ofstream ferr("abm_error_analysis.txt");
	ferr << "N0\tNf\tdt\tB\tEb\ttsys\n";
	
	for (int i=3; i<12; ++i){
		int N0 = pow(2,i);
		cout << "N0 = " << N0 << endl;
		
		Species<Daphnia> spp;
		Environment E;

		Solver S(SOLVER_ABM);
		S.control.abm_n0 = N0;

		S.setEnvironment(&E);
		S.addSpecies({100}, {0}, {1}, {false}, &spp, 0, -1);
		S.addSystemVariables({E.K});  // this can be done either before or after addSpecies()

		S.initialize();
		//S.print();
		
		auto t1 = high_resolution_clock::now();
		for (double t=0.05; t <= 100; t=t+0.5) {
			S.step_to(t);
		}
	    auto t2 = high_resolution_clock::now();
		duration<double, std::milli> ms_double = t2 - t1;
	
		
		double B = S.state_integral([&S](int i, double t){return S.species_vec[0]->getX(i)[0];}, S.current_time, 0);
		ferr << N0 << "\t" << S.species_vec[0]->xsize() << "\t" << 0 << "\t" << B << "\t" << fabs(B-1.298077)/1.298077 << "\t" << ms_double.count() << "\n";
	}

	ferr.close();
}	


int main(){
	run_fmu();
	run_ifmu();
	run_ebt();
	run_iebt();
	run_cm();
	run_icm();
	run_abm();

	return 0;
}


