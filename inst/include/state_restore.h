#ifndef PLANT_FATE_IO_STATE_RESTORE_H_
#define PLANT_FATE_IO_STATE_RESTORE_H_

#include <fstream>
#include <string>

#include <utils/initializer.h>
#include "trait_evolution.h"
#include "pspm_interface.h"

void saveState(Solver * S, std::string configfilename);
void restoreState(Solver * S, std::string configfilename);


#endif

