#include <iostream>
#include <iomanip>
#include <fstream>
#include "phydro.h"
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
	//double co2  = 400;          // ppm
	double elv  = 0;            // m.a.s.l.
	double fapar = 0.7;         // fractioni
	double rdark = 0.00;
	double psi_soil = 0;
	double pa = phydro::calc_patm(elv);

	phydro::ParCost par_cost(0.1, 1);
	phydro::ParPlant par_plant(3e-17, -2, 2);

	cout << setw(10) << "co2    " << "\t";
	cout << setw(10) << "jmax   " << "\t";
	cout << setw(10) << "dpsi   " << "\t";
	cout << setw(10) << "gs     " << "\t";
	cout << setw(10) << "a      " << "\t";
	cout << setw(10) << "ci     " << "\t";
	cout << setw(10) << "chi    " << "\t";
	cout << setw(10) << "vcmax  " << "\n";

	ifstream fin("tests/test_data/co2.tsv");
	double nerr = 0; int count = 0; 

	for (auto co2 : seq(150, 800, 50)){

		auto res = phydro::phydro_analytical(tc, tc, ppfd, ppfd/2, vpd, co2, pa, fapar, kphio, psi_soil, rdark, 3.0, par_plant, par_cost);
		
		cout << setw(10) <<  co2       << "\t"; cout.flush();
		cout << setw(10) <<  res.jmax  << "\t";
		cout << setw(10) <<  res.dpsi  << "\t";
		cout << setw(10) <<  res.gs    << "\t";
		cout << setw(10) <<  res.a     << "\t";
		cout << setw(10) <<  res.ci    << "\t";
		cout << setw(10) <<  res.chi   << "\t";
		cout << setw(10) <<  res.vcmax << "\n"; cout.flush();
	
		double d; 
		fin >> d; // first fin is psis
		fin >> d; nerr += err(res.jmax, d);  ++count;
		fin >> d; nerr += err(res.dpsi, d);  ++count;
		fin >> d; nerr += err(res.gs, d);    ++count;
		fin >> d; nerr += err(res.a, d);     ++count;
		fin >> d; nerr += err(res.ci, d);    ++count;
		fin >> d; nerr += err(res.chi, d);   ++count;
		fin >> d; nerr += err(res.vcmax, d); ++count;
		
	}
	
	cout << "Average error = " << nerr/count << "\n";
	
	if (abs (nerr/count) < 5e-5) return 0;
	else return 1;
}


