#ifndef PLANT_FATE_IO_STATE_RESTORE_H_
#define PLANT_FATE_IO_STATE_RESTORE_H_

#include <fstream>
#include <string>

#include <utils/initializer_v2.h>
#include "adaptive_species.h"
#include "pspm_interface.h"

void saveState(Solver * S, std::string state_outfile, std::string config_outfile, std::string params_file);
void restoreState(Solver * S, std::string state_infile, std::string config_infile);


#endif

