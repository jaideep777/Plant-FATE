#include <iomanip>
#include <fstream>
#include <random>

#include "life_history.h"
using namespace std;

int main(){

	cout << setprecision(12);
	
	plant::Plant P;
	P.initFromFile("tests/params/p_test_v2.ini");

	P.traits.print();

	cout << P.get_evolvableTraits({"lma", "wood_density", "hmat", "p50_xylem"}) << '\n';
	P.set_evolvableTraits({"lma", "wood_density", "hmat", "p50_xylem"}, {0.18, 800, 35, -4});

	P.traits.print();

	return 0;
}
