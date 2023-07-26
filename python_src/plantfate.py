from plantFATE import Simulator as sim

import pandas as pd

import numpy as np
from datetime import datetime

class PlantFATERunner: 

    iniFile = None
    plantFATE_model = None

    saveOutputs = True

    environment = pd.DataFrame(columns=['date', 'tair', 'ppfd_max', 'ppfd', 'vpd', 'elv', 'co2', 'swp'])
    emergentProps = pd.DataFrame()
    speciesProps = pd.DataFrame

    def __init__(self, param_file):
        self.plantFATE_model = sim("params/p_daily.ini")
        # self.soil_file = pd.read_csv(soil_file)
        

    def init(self, tstart, tend):
        self.plantFATE_model.init(tstart, tend)


    def init(self, tstart, temp0, soil_water_potentials0, vpd_0, par_0):
        photosynthetically_active_radiation = par_0 * 2.15
        newclim = self.plantFATE_model.Clim(temp0,
                                                photosynthetically_active_radiation * 4,
                                                photosynthetically_active_radiation,
                                                vpd_0 * 1000,
                                                np.nan,
                                                np.nan,
                                                soil_water_potentials0[0])
        
        datestart = datetime.strptime(tstart, "%Y-%m-%d")
        datediff = datetime.strptime(tstart, "%Y-%m-%d") - datetime(datestart.year, 1, 1)
        tstart = datestart.year + datediff.days/365
        self.plantFATE_model.init(tstart, newclim)


    def runstep(self, soil_water_potentials, vapour_pressure_deficit, photosynthetically_active_radiation,
                      temperature):
        photosynthetically_active_radiation = photosynthetically_active_radiation * 2.15
        self.plantFATE_model.update_environment(temperature,
                                                photosynthetically_active_radiation * 4,
                                                photosynthetically_active_radiation,
                                                vapour_pressure_deficit * 1000,
                                                np.nan,
                                                np.nan,
                                                soil_water_potentials[0])
        self.plantFATE_model.simulate_step()

        self.saveEnvironment(self)
        self.saveEmergentProps(self)

        trans = self.plantFATE_model.props.trans
        # return evapotranspiration, soil_specific_depletion_1, soil_specific_depletion_2, soil_specific_depletion_3
        return trans, 0, 0, 0

    def saveEnvironment(self):
        e = pd.DataFrame({'date' : [self.plantFATE_model.E.tcurrent],
                          'tair' : [self.plantFATE_model.E.weightedAveClim.tc], 
                          'ppfd_max' : [self.plantFATE_model.E.weightedAveClim.ppfd_max], 
                          'ppfd' : [self.plantFATE_model.E.weightedAveClim.ppfd], 
                          'vpd' : [self.plantFATE_model.E.weightedAveClim.vpd] ,
                          'elv' : [self.plantFATE_model.E.weightedAveClim.elv],
                          'co2' : [self.plantFATE_model.E.weightedAveClim.co2], 
                          'swp' : [self.plantFATE_model.E.weightedAveClim.swp]})
        
        self.environment.append(e)

    def saveEmergentProps(self):
        e = pd.DataFrame({'date' : [self.plantFATE_model.tcurrent],
                          'trans' : [self.plantFATE_model.props.trans], 
                          'gs' : [self.plantFATE_model.props.gs], 
                          'gpp' : [self.plantFATE_model.props.gpp], 
                          'lai' : [self.plantFATE_model.props.lai] ,
                          'npp' : [self.plantFATE_model.props.npp],
                          'cc_est' : [self.plantFATE_model.props.cc_est], 
                          'croot_mass' : [self.plantFATE_model.props.croot_mass],
                          'froot_mass' : [self.plantFATE_model.props.froot_mass],
                          'lai_vert' : [self.plantFATE_model.props.lai_vert],
                          'leaf_mass' : [self.plantFATE_model.props.leaf_mass],
                          'resp_auto' : [self.plantFATE_model.props.resp_auto],
                          'stem_mass' : [self.plantFATE_model.props.stem_mass]})
        
        self.emergentProps.append(e)
    