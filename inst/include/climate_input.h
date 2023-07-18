#ifndef PLANT_FATE_ENV_CLIMATE_INPUT_H_
#define PLANT_FATE_ENV_CLIMATE_INPUT_H_

#include "climate.h"


namespace env{

class ClimateInput{
    public:

    Clim currentClim; // return the 
    Clim weightedAveClim;

    double tcurrent;
    double ave_window = 14;

    ClimateInput();
    ClimateInput(Clim &climObj, double t0, double _ave_window);

    void updateEnvironment();
    void updateClim(Clim &newClim, double tnew);
    void print_line();

    private:

    double movingAverageRungeKutta4(double xt1, double xt2, double h, double yt1);
    void checkIfNaN(Clim &newClimObj);
    Clim computeNewAverage(Clim &newClim, double tnew);
};

}


#endif
