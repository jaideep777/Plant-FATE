#include <iostream>
#include <iomanip>
#include <fstream>
#include "phydro.h"
#include <chrono>
using namespace std;


double err(double x, double ref){
	return std::min(abs(x-ref), abs(x/ref-1));
}

int check(double x, double ref, double err=1e-4){
	//cout << "err: " << abs(x-ref) << " " << abs(x/ref-1)<< "\n";
	//cout << "comp: " << x << " " << ref << "|" << abs(x-ref) << "\n";
	if (abs(x-ref) < err || abs(x/ref-1) < err) return 0;
	else return 1;
}


vector<double> seq(double start, double end, int length){
	vector<double> x;
	for (int i=0; i<length; ++i) x.push_back(start + double(i)*(end-start)/(length-1));
	return x;
}

vector<double> lseq(double start, double end, int length){
	vector<double> x;
	for (int i=0; i<length; ++i) x.push_back(exp(log(start) + double(i)*(log(end)-log(start))/(length-1)));
	return x;
}

int main(){

	double kphio = 0.087;        // quantum yield efficiency
	// c_molmass <- 12.0107 # molar mass, g / mol

	// Define environmental conditions
	double tc = 25;             // temperature, deg C
	double ppfd = 300;          // umol/m2/s
	double vpd  = 1000;         // Pa
	double co2  = 400;          // ppm
	double elv  = 0;            // m.a.s.l.
	double fapar = 0.7;         // fractioni
	double rdark = 0.02;
	double pa = phydro::calc_patm(elv);
	
	phydro::ParCost par_cost(0.1, 1);
	phydro::ParPlant par_plant(3e-17, -2, 2);
	phydro::ParControl par_control;

	double nerr = 0; int count = 0;

	double time_num = 0, time_ana = 0, time_num_apx = 0, time_ana_apx = 0, time_num_apx2 = 0, time_ana_apx2 = 0, time_num_qng = 0, time_ana_qng = 0;
	for (auto psi_soil : seq(-6, 0, 1000)){

		par_control.gs_method = phydro::GS_IGF;
		
		auto t1 = std::chrono::high_resolution_clock::now();
		phydro::phydro_numerical(tc, tc, ppfd, ppfd/2, vpd, co2, pa, fapar, kphio, psi_soil, rdark, 3.0, par_plant, par_cost, par_control);
		auto t2 = std::chrono::high_resolution_clock::now();
		time_num += (std::chrono::duration<double, std::micro> (t2 - t1)).count();
			
		auto t3 = std::chrono::high_resolution_clock::now();
		phydro::phydro_analytical(tc, tc, ppfd, ppfd/2, vpd, co2, pa, fapar, kphio, psi_soil, rdark, 3.0, par_plant, par_cost, par_control);
		auto t4 = std::chrono::high_resolution_clock::now();
		time_ana += (std::chrono::duration<double, std::micro> (t4 - t3)).count();
		
		
		par_control.gs_method = phydro::GS_QNG;
		
		auto t9 = std::chrono::high_resolution_clock::now();
		phydro::phydro_numerical(tc, tc, ppfd, ppfd/2, vpd, co2, pa, fapar, kphio, psi_soil, rdark, 3.0, par_plant, par_cost, par_control);
		auto t10 = std::chrono::high_resolution_clock::now();
		time_num_qng += (std::chrono::duration<double, std::micro> (t10 - t9)).count();
			
		auto t11 = std::chrono::high_resolution_clock::now();
		phydro::phydro_analytical(tc, tc, ppfd, ppfd/2, vpd, co2, pa, fapar, kphio, psi_soil, rdark, 3.0, par_plant, par_cost, par_control);
		auto t12 = std::chrono::high_resolution_clock::now();
		time_ana_qng += (std::chrono::duration<double, std::micro> (t12 - t11)).count();
		
		
		par_control.gs_method = phydro::GS_APX;

		auto t5 = std::chrono::high_resolution_clock::now();
		phydro::phydro_numerical(tc, tc, ppfd, ppfd/2, vpd, co2, pa, fapar, kphio, psi_soil, rdark, 3.0, par_plant, par_cost, par_control);
		auto t6 = std::chrono::high_resolution_clock::now();
		time_num_apx += (std::chrono::duration<double, std::micro> (t6 - t5)).count();
			
		auto t7 = std::chrono::high_resolution_clock::now();
		phydro::phydro_analytical(tc, tc, ppfd, ppfd/2, vpd, co2, pa, fapar, kphio, psi_soil, rdark, 3.0, par_plant, par_cost, par_control);
		auto t8 = std::chrono::high_resolution_clock::now();
		time_ana_apx += (std::chrono::duration<double, std::micro> (t8 - t7)).count();

		par_control.gs_method = phydro::GS_APX2;
		
		auto t13 = std::chrono::high_resolution_clock::now();
		phydro::phydro_numerical(tc, tc, ppfd, ppfd/2, vpd, co2, pa, fapar, kphio, psi_soil, rdark, 3.0, par_plant, par_cost, par_control);
		auto t14 = std::chrono::high_resolution_clock::now();
		time_num_apx2 += (std::chrono::duration<double, std::micro> (t14- t13)).count();
			
		auto t15 = std::chrono::high_resolution_clock::now();
		phydro::phydro_analytical(tc, tc, ppfd, ppfd/2, vpd, co2, pa, fapar, kphio, psi_soil, rdark, 3.0, par_plant, par_cost, par_control);
		auto t16 = std::chrono::high_resolution_clock::now();
		time_ana_apx2 += (std::chrono::duration<double, std::micro> (t16 - t15)).count();

	}
	
	cout << "                | " << "Analytical  \t Numerical"   << "\n";
	cout << "----------------+---------------------------------------\n";
	cout << "Avg time (IGF)  | " << time_ana/1000     << " us\t" << time_num/1000      << " us\n";
	cout << "Avg time (QNG)  | " << time_ana_qng/1000 << " us\t" << time_num_qng/1000  << " us\n";
	cout << "Avg time (APX)  | " << time_ana_apx/1000 << " us\t" << time_num_apx/1000  << " us\n";
	cout << "Avg time (APX2) | " << time_ana_apx2/1000 << " us\t" << time_num_apx2/1000 << " us\n";
	cout << "-------------------------------------------------------\n";
	cout << "Speedup (IGF)  = " << time_num/time_ana          << " x\n";
	cout << "Speedup (QNG)  = " << time_num_qng/time_ana_qng  << " x\n";
	cout << "Speedup (APX)  = " << time_num_apx/time_ana_apx  << " x\n";
	cout << "Speedup (APX2) = " << time_num_apx2/time_ana_apx2 << " x\n";
	
	return 0;
}


