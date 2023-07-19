#include <math.h>

#include "climate_input.h"

namespace env{

ClimateInput::ClimateInput(){
    weightedAveClim = Clim();
    currentClim = Clim();
    tcurrent = 0;
}

ClimateInput::ClimateInput(Clim &climObj, double t0, double _ave_window){
    weightedAveClim = climObj;
    currentClim = climObj;
    tcurrent = t0;
    ave_window = _ave_window;
}

void ClimateInput::updateEnvironment(){}

void ClimateInput::updateClim(Clim &newClim, double tnew){
    checkIfNaN(newClim);

    // weightedAveClim = computeNewAverage(newClim, tnew);
    weightedAveClim = newClim;
    currentClim = newClim;
    tcurrent = tnew;
}

Clim ClimateInput::computeNewAverage(Clim &newClim, double tnew){
    double h = tnew - tcurrent;

    Clim newWeightedAve;
    
    newWeightedAve.tc = movingAverageRungeKutta4(currentClim.tc, newClim.tc, h, weightedAveClim.tc);
    newWeightedAve.ppfd_max = movingAverageRungeKutta4(currentClim.ppfd_max, newClim.ppfd_max, h, weightedAveClim.ppfd_max);
    newWeightedAve.ppfd = movingAverageRungeKutta4(currentClim.ppfd, newClim.ppfd, h, weightedAveClim.ppfd);
    newWeightedAve.vpd = movingAverageRungeKutta4(currentClim.vpd, newClim.vpd, h, weightedAveClim.vpd);
    newWeightedAve.co2 = movingAverageRungeKutta4(currentClim.co2, newClim.co2, h, weightedAveClim.co2);
    newWeightedAve.elv = movingAverageRungeKutta4(currentClim.elv, newClim.elv, h, weightedAveClim.elv);
    newWeightedAve.swp = movingAverageRungeKutta4(currentClim.swp, newClim.swp, h, weightedAveClim.swp);

    return newWeightedAve;
}

double ClimateInput::movingAverageRungeKutta4(double xt1, double xt2, double h, double yt1){
    double h2 = h/2;

    double xt15 = (xt1 + xt2)/2; //linear interpolation

    double k1 = 1/ave_window * xt1 - 1/ave_window * yt1;
    double k2 = 1/ave_window * xt15 - 1/ave_window * (yt1 + h2 * k1);
    double k3 = 1/ave_window * xt15 - 1/ave_window * (yt1 + h2 * k2);
    double k4 = 1/ave_window * xt2 - 1/ave_window * (yt1 + h * k3);

    double yt2 = yt1 + h/6 * (k1 + 2 * k2 + 2 * k3 + k4);
    yt2 = xt2;

    return yt2;
}

void ClimateInput::checkIfNaN(Clim &newClimObj){
    if(isnan(newClimObj.tc)){
        newClimObj.tc = currentClim.tc;
    }
    if(isnan(newClimObj.ppfd_max)){
        newClimObj.ppfd_max = currentClim.ppfd_max;
    }
    if(isnan(newClimObj.ppfd)){
        newClimObj.ppfd = currentClim.ppfd;
    }
    if(isnan(newClimObj.vpd)){
        newClimObj.vpd = currentClim.vpd;
    }
    if(isnan(newClimObj.co2)){
        newClimObj.co2 = currentClim.co2;
    }
    if(isnan(newClimObj.elv)){
        newClimObj.elv = currentClim.elv;
    }
    if(isnan(newClimObj.swp)){
        newClimObj.swp = currentClim.swp;
    }
}

void ClimateInput::print_line(){
    std::cout.flush();
    std::cout << "Climate at t = " << tcurrent;
    std::cout.flush();
	std::cout << " | T/Im/I/D/CO2/Z/Ps = " << currentClim.tc << " " << currentClim.ppfd_max << " " << currentClim.ppfd << " " << currentClim.vpd << " " << currentClim.co2 << " " << currentClim.elv << " " << currentClim.swp << "\n"; 
    std::cout.flush();
    std::cout << "Weighted Climate at t = " << tcurrent;
    std::cout.flush();
	std::cout << " | T/Im/I/D/CO2/Z/Ps = " << weightedAveClim.tc << " " << weightedAveClim.ppfd_max << " " << weightedAveClim.ppfd << " " << weightedAveClim.vpd << " " << weightedAveClim.co2 << " " << weightedAveClim.elv << " " << weightedAveClim.swp << "\n"; 
    std::cout << std::endl;
    std::cout.flush();
}

}
