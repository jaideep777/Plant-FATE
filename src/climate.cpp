#include <iomanip>
#include <cmath>
#include <stdexcept>

#include "climate.h"
using namespace std;

namespace env{

void Climate::set_elevation(double _elv){
	clim_inst.set_elevation(_elv);
	clim_acclim.set_elevation(_elv);
}

void Climate::set_co2(double _co2){
	clim_inst.co2 = _co2;
	clim_acclim.co2 = _co2;
}

void Climate::print(double t){
	cout << "Current climate (inst / acclim):\n";
	cout << "   tc      = " << clim_inst.tc       << " / " << clim_acclim.tc       << '\n';
	cout << "   ppfd    = " << clim_inst.ppfd     << " / " << clim_acclim.ppfd     << '\n';
	cout << "   vpd     = " << clim_inst.vpd      << " / " << clim_acclim.vpd      << '\n';
	cout << "   co2     = " << clim_inst.co2      << " / " << clim_acclim.co2      << '\n';
	cout << "   elv     = " << clim_inst.elv      << " / " << clim_acclim.elv      << '\n';
	cout << "   swp     = " << clim_inst.swp      << " / " << clim_acclim.swp      << '\n';
}

void Climate::print_line(double t){
	cout << "Current climate (t=" << t << "): ";
	cout << clim_inst.tc  << " ";
	cout << clim_inst.ppfd << " ";
	cout << clim_inst.vpd << " ";
	cout << clim_inst.co2 << " ";
	cout << clim_inst.elv << " ";
	cout << clim_inst.swp << " ";
	cout << '\n';
}

} // namespace env


