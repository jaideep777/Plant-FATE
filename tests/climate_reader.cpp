#include <iostream>
#include <fstream>

#include "climate.h"

using namespace std;

int main(){
	
	env::Climate C;

	C.metFile = "tests/data/MetData_AmzFACE_Monthly_2000_2015_PlantFATE.csv";
	C.co2File = "tests/data/CO2_ELE_AmzFACE2000_2100.csv";

	C.init();
	C.print(0);

	ofstream fout("climate.txt");
	for (double t = 2015; t < 2017; t += 1/24.0){
		C.updateClimate(t-2000);
		C.print(t-2000);
		fout << t << "\t" << C.clim.tc << "\t" << C.clim.vpd << "\t" << C.clim.ppfd << "\t" << C.clim.swp << "\n";
	}
	fout.close();

	return 0;

}


//dat = read.delim("/home/jaideep/codes/plant_fate_ppa/climate.txt", header=F)
//dato = read.csv("/home/jaideep/codes/plant_fate_ppa/tests/data/MetData_AmzFACE_Monthly_2000_2015_PlantFATE.csv", header=T)
//dato$t = dato$Year + (dato$Month-1)/12

//par(mfrow=c(2,2), mar=c(4,4,1,1), oma=c(1,1,1,1))
//plot(dat$V2~dat$V1, type="l", col="black")
//points(dato$Temp~dato$t, pch=20, type="p", col="red")
//plot(dat$V3~dat$V1, type="l", col="black")
//points(dato$VPD~dato$t, pch=20, type="p", col="red")
//plot(dat$V4~dat$V1, type="l", col="black")
//points(dato$PAR~dato$t, pch=20, type="p", col="red")
//plot(dat$V5~dat$V1, type="l", col="black")
//points(dato$SWP~dato$t, pch=20, type="p", col="red")

