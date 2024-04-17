import os.path
import sys
import EZPaths
import classRunComputing
import runDataClass

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


# Test Run
# current_run = "testRun1"
# current_bp = astroModels.bp_testRun1
# current_var_params = ["p_temp","p_mag"]
# current_geo_grid_names = ["ModelA", "ModelB"]
# current_geo_grid_values = [(["a"], [str(.3)]),(["a"], [str(.9)])]
# action = {
#     "var": "nu0",
#     "start": 3.00e+10,
#     "stop": 6.00e+10,
#     "step": 1.00e+10,
#     "images": True
# }

# Test Run 2
# current_run = "testRun2"
# current_bp = astroModels.bp_testRun1
# current_var_params = ["p_temp","p_mag"]
# current_geo_grid_names = ["ModelA", "ModelB"]
# current_geo_grid_values = [(["a"], [str(.3)]),(["a"], [str(.9)])]
# action = {
#     "var": "nu0",
#     "start": 3.00e+10,
#     "stop": 6.00e+10,
#     "step": 1.00e+10,
#     "images": True
# }
# isNormalized = True
# blurr_policy = True
# action = {
#     "var": "nu0",
#     "start": 670e9,
#     "stop": 700e9,
#     "step": 20e9,
#     "images": True
# }
# isNormalized = True
# blurr_policy = True

# Solo Run
# current_run = "soloRun1"
# current_bp = astroModels.bp_soloRun1
# current_var_params = ["p_mag"]
# current_geo_grid_names = ["ModelA","ModelB"]
# current_geo_grid_values = [(["a"], [str(.3)]),(["a"], [str(.9)])]
# action = {
#     "var": "nu0",
#     "start": 10e9,
#     "stop": 700e9,
#     "step": 20e9,
#     "images": True
# }

# Solo Run 2
# current_run = "soloRun2"
# current_bp = astroModels.bp_soloRun2
# current_var_params = []
# current_geo_grid_names = ["ModelB"]
# current_geo_grid_values = [(["a"], [str(.9)])]
# action = {
#     "var": "nu0",
#     "start": 670e9,
#     "stop": 700e9,
#     "step": 20e9,
#     "images": True
# }
# isNormalized = True
# blurr_policy = True


#Full Run
current_run = "run1"
current_bp = astroModels.bp_run1
current_var_params = ["p_temp", "p_mag"]
current_geo_grid_names = ["ModelA", "ModelB"]
current_geo_grid_values = [(["a"], [str(.3)]),(["a"], [str(.9)])]
action = {
    "var": "nu0",
    "start": 10e9,
    "stop": 700e9,
    "step": 20e9,
    "images": True
}
#     "var": "nu0",
#     "start": 670e9,
#     "stop": 700e9,
#     "step": 20e9,
#     "images": True
# }
isNormalized = True
blurr_policy = True

#
# sub_paths, all_intensity_models, total_models_count, run_type, variable_param_ranges, constant_params =\
#     fileloading.runsInit(current_run,current_bp,current_var_params)
# bigRunComputing.createGeoGrid(sub_paths, current_geo_grid_names, current_run)
# bigRunComputing.creatIntensityGrid(sub_paths, current_run, current_geo_grid_names, current_geo_grid_values,
#                                    all_intensity_models, current_var_params, constant_params,
#                                    total_models_count,action)
# bigRunComputing.graphCreation(sub_paths,current_run,action,2)
# bigRunComputing.surfacePlot(sub_paths,current_bp,action,current_var_params,current_geo_grid_names)


# bigRun = classRunComputing.BigRuns(current_run,current_bp,current_var_params,
#                                    current_geo_grid_values,current_geo_grid_names)
# # bigRun.createGeoGrid()
# bigRun.creatIntensityGrid(action,isNormalized,blurr_policy)
#
# bigRun.blurrGraphCreation(action)
# bigRun.graphCreation(action)
#

run = runDataClass.testRun2

bigRun = classRunComputing.BigRuns(
    run.getRunName(),
    run.getBrightparams(),
    run.getBPVarNames(),
    run.getGeoGrid(),
    run.getGeoGridNames()
)
""" Geo Model_________________________________________________________"""

# bigRun.createGeoGrid()

""" Intensity Grid_________________________________________________________"""

# bigRun.creatIntensityGrid(
#     run.getAction(),
#     run.getIsNormalized(),
#     run.getBlurrPolicy()
# )
""" Blurr Graph Creation_________________________________________________________"""

bigRun.blurrGraphCreation(
    run.getAction()
)
""" _________________________________________________________"""

bigRun.graphCreation(
    run.getAction()
)
""" _________________________________________________________"""


