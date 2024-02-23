#include <iomanip>
#include <cmath>
#include <stdexcept>

#include "climate.h"
using namespace std;

namespace env{

void Climate::set_elevation(double _elv){
	clim.set_elevation(_elv);
	clim_acclim.set_elevation(_elv);
}

void Climate::print(double t){
	cout << "Current climate (inst / acclim):\n";
	cout << "   tc      = " << clim.tc       << " / " << clim_acclim.tc       << '\n';
	cout << "   ppfdmax = " << clim.ppfd_max << " / " << clim_acclim.ppfd_max << '\n';
	cout << "   ppfd    = " << clim.ppfd     << " / " << clim_acclim.ppfd     << '\n';
	cout << "   vpd     = " << clim.vpd      << " / " << clim_acclim.vpd      << '\n';
	cout << "   co2     = " << clim.co2      << " / " << clim_acclim.co2      << '\n';
	cout << "   elv     = " << clim.elv      << " / " << clim_acclim.elv      << '\n';
	cout << "   swp     = " << clim.swp      << " / " << clim_acclim.swp      << '\n';
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


