from aart_func import *
import params
import bareFileLoading as bfileloading
import fileloading
from params import * # The file params.py contains all the relevant parameters for the simulations
import importlib

print(spin_case)
current_geo_model = fileloading.totalModelNametoGridModel("ModelA222",3)
bfileloading.loadGeoModel(current_geo_model, "PRUN")()
print(spin_case)