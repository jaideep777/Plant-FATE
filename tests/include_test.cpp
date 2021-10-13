#include "plant_params.h"
#include "plant_geometry.h"
#include "initializer.h"
using namespace std;

int main(){

	plant::PlantParameters par;

	par.initFromFile("tests/params/p.ini");
	par.print();

	cout << "***************************************\n";

	io::Initializer I("tests/params/p.ini");
	I.readFile();
	I.print();
	double a = I.get<double>("testVal"); 
	cout << "number in string testVal = " << a << "\n";

	return 0;
}

