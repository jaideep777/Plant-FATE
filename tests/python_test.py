import numpy as np
import pandas as pd
import psutil

from pypfate import Patch as patch
from pypfate import Clim

## Read in dataframe with climate data
climate_data = pd.read_csv("~/drought-MIP/test_data.csv")
print(climate_data.head())
print(climate_data.loc[0,"Time"])

## Initialise plantfate patch
PFpatch = patch(str("tests/params/python_p.ini"))
PFpatch.init(climate_data.loc[0,"Time"], climate_data.loc[0,"Time"] + 1000)

PFpatch.update_climate(368.9,
                              climate_data.loc[0,"temp"],
                              climate_data.loc[0,"VPD"],
                              climate_data.loc[0,"PPFD"],
                              climate_data.loc[0,"SWP"],
                              climate_data.loc[0,"NR"])

## Run loop with plantfate - save yearly memory details
# for i in climate_data.loc[1:(len(climate_data) - 1),"Time"] :
mem = []
#
for i in range(1,len(climate_data)):
# for i in range(1,366):
    PFpatch.simulate_to(climate_data.loc[i,"Time"])
    PFpatch.update_climate(368.9,
                              climate_data.loc[i,"temp"],
                              climate_data.loc[i,"VPD"],
                              climate_data.loc[i,"PPFD"],
                              climate_data.loc[i,"SWP"],
                              climate_data.loc[i,"NR"])
    process = psutil.Process()
    mem.append(float(process.memory_info().rss) / (1024 ** 3))

df = pd.DataFrame(data={"Memory": mem})
df.to_csv("memory_use.csv")
#
## close plantfate
PFpatch.close()