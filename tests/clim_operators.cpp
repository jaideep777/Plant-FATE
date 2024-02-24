#include <iostream>
#include <fstream>

#include "climate.h"

using namespace std;

int main(){
	
	env::Clim C;
	C.print();

	C *= 0.1;
	C.print();

	C *= 10;

	env::Clim C2 = C*0.1;
	C2.print();

	env::Clim Cnext = C;
	Cnext *= 1.5;
	Cnext.print();

	env::Clim C3 = C;
	C3 += (1-exp(-0.1))*(Cnext - C3);
	C3.print();

	double tc_expected = C.tc + (1-exp(-0.1))*(C.tc*1.5 - C.tc);
	double ppfd_expected = C.ppfd + (1-exp(-0.1))*(C.ppfd*1.5 - C.ppfd);
	double rn_expected = C.rn + (1-exp(-0.1))*(C.rn*1.5 - C.rn);
	double vpd_expected = C.vpd + (1-exp(-0.1))*(C.vpd*1.5 - C.vpd);
	double co2_expected = C.co2 + (1-exp(-0.1))*(C.co2*1.5 - C.co2);
	double elv_expected = C.elv + (1-exp(-0.1))*(C.elv*1.5 - C.elv);
	double swp_expected = C.swp + (1-exp(-0.1))*(C.swp*1.5 - C.swp);
	double vwind_expected = C.vwind + (1-exp(-0.1))*(C.vwind*1.5 - C.vwind);
	double pa_expected = C.pa + (1-exp(-0.1))*(C.pa*1.5 - C.pa);

	if (fabs(tc_expected - C3.tc) > 1e-6) return 1;
	if (fabs(ppfd_expected - C3.ppfd) > 1e-6) return 1;
	if (fabs(rn_expected - C3.rn) > 1e-6) return 1;
	if (fabs(vpd_expected - C3.vpd) > 1e-6) return 1;
	if (fabs(co2_expected - C3.co2) > 1e-6) return 1;
	if (fabs(elv_expected - C3.elv) > 1e-6) return 1;
	if (fabs(swp_expected - C3.swp) > 1e-6) return 1;
	if (fabs(vwind_expected - C3.vwind) > 1e-6) return 1;
	if (fabs(pa_expected - C3.pa) > 1e-6) return 1;

	return 0;

}

