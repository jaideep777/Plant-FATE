#include <iostream>
#include <iomanip>
#include <fstream>
#include "phydro.h"
using namespace std;


double err(double x, double ref){
	return std::min(abs(x-ref), abs(x/ref-1));
}

int check(double x, double ref, double err=1e-6){
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
	phydro::ParEnv my_env(25.0, 101325.0, 2000.0, 300.0, 5.0);
	phydro::ParEnv my_env2(20.0, 101325.0, 1000.0, 200.0);

	my_env.print();
	my_env2.print();
	return 0;

}

