#include <iomanip>
#include "exp_average.h"
#include <fstream>

using namespace std;

int main(){
	double days_per_year = 365.2524;
	ExpAverager<double> M(2000,25, 14/days_per_year);

	ofstream fout("exp_avg.txt");	
	fout << "t\tf\tf_avg\n";
	for (int i=0; i<1000; ++i){
		double t = 2000 + i/days_per_year;
		double f = 25 + sin(i*2*M_PI/days_per_year) + 0.3*(1-2*double(rand())/RAND_MAX);
		M.update(t, f);
		fout << setprecision(12)
		     << t << "\t"
		     << f << "\t"
			 << M.get() << "\n";
	}

	return 0;
}

