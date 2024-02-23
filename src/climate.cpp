#include <iomanip>
#include <cmath>
#include <stdexcept>

#include "climate.h"
using namespace std;

namespace env{

void Climate::set_elevation(double _elv){
	clim.set_elevation(_elv);
}

void Climate::print(double t){
	cout << "Current climate:\n";
	cout << "   tc      = " << clim.tc  << '\n';
	cout << "   ppfdmax = " << clim.ppfd_max << '\n';
	cout << "   ppfd    = " << clim.ppfd << '\n';
	cout << "   vpd     = " << clim.vpd << '\n';
	cout << "   co2     = " << clim.co2 << '\n';
	cout << "   elv     = " << clim.elv << '\n';
	cout << "   swp     = " << clim.swp << '\n';
}

void Climate::print_line(double t){
	cout << "Current climate (t=" << t << "): ";
	cout << clim.tc  << " ";
	cout << clim.ppfd_max << " ";
	cout << clim.ppfd << " ";
	cout << clim.vpd << " ";
	cout << clim.co2 << " ";
	cout << clim.elv << " ";
	cout << clim.swp << " ";
	cout << '\n';
}

} // namespace env


