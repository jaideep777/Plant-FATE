#include <iostream>
#include <fstream>

#include "climate.h"

using namespace std;

int main(){
	
	env::Climate C;

	C.metFile = "tests/data/MetData_AmzFACE_Monthly_2000_2001_PlantFATE.csv";
	C.co2File = "tests/data/CO2_AMB_AmzFACE2000_2100.csv";
	
	C.update_met = true;
	C.update_co2 = true;

	C.init();
	C.print(0);

	ofstream fout("climate.txt");
	C.t0 = 1950;
	for (double t = 1950; t < 2025; t += 1/120.0){
		C.updateClimate(t);
		C.print(t);
		fout << t << "\t" << C.clim.tc << "\t" << C.clim.vpd << "\t" << C.clim.ppfd << "\t" << C.clim.swp << "\t" << C.clim.co2 << "\n";
	}
	fout.close();

	return 0;

}

//dat = read.delim("~/codes/tmodel_cpp/climate.txt", header=F)
//dato = read.csv("~/codes/tmodel_cpp/tests/data/MetData_AmzFACE_Monthly_2000_2001_PlantFATE.csv", header=T)
//dato$t = dato$Year + (dato$Month-1)/12
//dato = rbind(dato, dato, dato, dato, dato, dato)
//dato$Year = c(rep(2000,12), rep(2001,12), 
//              rep(2002,12), rep(2003,12),
//              rep(2004,12), rep(2005,12),
//              rep(2006,12), rep(2007,12),
//              rep(2008,12), rep(2009,12),
//              rep(2010,12), rep(2011,12))
//dato$t = dato$Year + (dato$Month-1)/12

//par(mfrow=c(4,1), mar=c(4,4,1,1), oma=c(1,1,1,1))
//plot(dat$V2~dat$V1, type="l", col="black")
//points(dato$Temp~dato$t, pch=20, type="p", col="red")
//plot(dat$V3~dat$V1, type="l", col="black")
//points(dato$VPD~dato$t, pch=20, type="p", col="red")
//plot(dat$V4~dat$V1, type="l", col="black")
//points(dato$PAR~dato$t, pch=20, type="p", col="red")
//plot(dat$V5~dat$V1, type="l", col="black")
//points(dato$SWP~dato$t, pch=20, type="p", col="red")
