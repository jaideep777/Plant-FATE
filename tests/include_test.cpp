#include "utils/rk4.h"
#include "utils/incbeta.h"
#include "plant_geometry.h"
#include "utils/initializer.h"
#include "utils/lambertw.h"
//#include "assimilation.h"
//
#include "trait_reader.h"

using namespace std;

int main(){
	
	plant::PlantParameters par;

	par.initFromFile("tests/params/p.ini");
	par.print();

	cout << "***************************************\n";

	io::Initializer I("tests/params/p.ini");
	I.readFile();
	I.print();
	double a = I.get<double>("outDir"); 
	cout << "number in string outDir = " << a << "\n";

	// lambert W test
	for (int i=0; i<21; ++i){
		double x = double(i);

		double w = lambertw0(x);
		if (fabs(w*exp(w) - x) > 1e-8) return 1;
	}

	// traits file test
	TraitsReader Tr;
	Tr.readFromFile("tests/data/trait1.csv");
	Tr.print();

	return 0;
}

