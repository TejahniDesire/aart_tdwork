import os.path
import sys
import EZPaths
sys.path.append(EZPaths.aartPath)
import bigRunComputing
import fileloading
import image_tools
import subprocess
import scipy
from matplotlib import ticker
from aart_func import *
from params import *  # The file params.py contains all the relevant parameters for the simulations
from astropy import units as u
import image_tools as tls
import numpy as np
import subprocess
import EZPaths
import os
from aart_func import *
from params import *
import importlib
import params
import astroModels
import fileloading


current_run = "testRun1"
current_bp = astroModels.bp_testRun1
current_var_params = ["p_temp","p_mag"]
current_geo_grid = ["ModelA", "ModelB"]
current_geo_grid_values = [(["a"], [str(.3)]),(["a"], [str(.9)])]
sub_paths, all_intent_models, model_name_string = fileloading.runsInit(current_run,current_bp,current_var_params)


# bigRunComputing.fileloading.loadGeoModel(current_geo_grid[0],current_run)
# bigRunComputing.createGeoGrid(sub_paths, current_geo_grid, current_run)
# bigRunComputing.creatIntensityGrid(sub_paths,current_geo_grid,current_run,all_intent_models,model_name_string,current_geo_grid_values)
bigRunComputing.graphCreation(sub_paths,current_run,bigRunComputing.kw_action,2)