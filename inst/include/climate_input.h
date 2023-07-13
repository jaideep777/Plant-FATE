#ifndef PLANT_FATE_ENV_CLIMATE_INPUT_H_
#define PLANT_FATE_ENV_CLIMATE_INPUT_H_

#include "climate.h"


namespace env{

class ClimateInput{
    public:

    Clim currentClim; // return the 
    Clim weightedAveClim;

    double tcurrent;
    double ave_window;



    void updateClim(Clim newClim, double tnew);


};


}


#endif
