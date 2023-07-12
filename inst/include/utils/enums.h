#ifndef ENUMS_H
#define ENUMS_H

#include <string>

using namespace std;

enum plant_solv_time_step {DAILY_STEP, YEARLY_STEP};

inline plant_solv_time_step timestep_map(std::string step){
    std::map<std::string, plant_solv_time_step> tmap = 
		{{"daily",  DAILY_STEP}, 
		 {"yearly",  YEARLY_STEP}};
    return tmap.at(step);
}

#endif