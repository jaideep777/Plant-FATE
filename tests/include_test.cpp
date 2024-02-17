#include "utils/rk4.h"
#include "utils/incbeta.h"
#include "plant_geometry.h"
#include "utils/initializer_v2.h"
#include "utils/lambertw.h"
//#include "assimilation.h"
//
#include "trait_reader.h"

using namespace std;

int main(){
	
	plant::PlantParameters par;

	par.initFromFile("tests/params/p_test_v2.ini");
	par.print();

	cout << "***************************************\n";

	io::Initializer I;
	I.parse("tests/params/p_test_v2.ini");
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
	Tr.readFromFile("tests/data/traits_example_canopy_understorey.csv");
	Tr.print();

	return 0;
}

