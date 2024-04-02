import os.path
import sys
import EZPaths
import classRunComputing
import smallRunComputing

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
#     "step": 1.00e+10
# }
#


# Full Run
# current_run = "run1"
# current_bp = astroModels.bp_run1
# current_var_params = ["p_temp", "p_mag"]
# current_geo_grid_names = ["ModelA", "ModelB"]
# current_geo_grid_values = [(["a"], [str(.3)]),(["a"], [str(.9)])]
# action = {
#     "var": "nu0",
#     "start": 10e9,
#     "stop": 700e9,
#     "step": 20e9
# }
# current_model = 'ModelB13'

# Solo Run
current_run = "soloRun1"
current_bp = astroModels.bp_soloRun1
current_var_params = []
current_geo_grid_names = ["ModelB"]
current_geo_grid_values = [(["a"], [str(.9)])]
action = {
    "var": "nu0",
    "start": 10e9,
    "stop": 700e9,
    "step": 20e9,
    "images": True
}
current_model = ['ModelB1']


save_paths = {
    'intVRad': '/scratch/gpfs/td6241/aart/bigRuns/' + current_run + '/Images/inensityVRadiiTEST/',
    'intVRad2': '/scratch/gpfs/td6241/aart/bigRuns/' + current_run + '/Images/inensityVRadiiTESTTWO/',
    'radVVarphi':'/scratch/gpfs/td6241/aart/bigRuns/' + current_run + '/Images/radVVarphi/'
}

for path in list(save_paths):
    fileloading.creatSubDirectory(save_paths[path])

'''_______________________________________________________________________'''
sub_paths, all_intensity_models, total_models_count, run_type, variable_param_ranges, constant_params =\
    fileloading.runsInit(current_run,current_bp,current_var_params)


smallRunComputing .playModel(sub_paths,save_paths, current_run,action, current_model, intent_grid_type=2)


bigRun = classRunComputing.BigRuns(current_run,current_bp,current_var_params,
                                   current_geo_grid_values,current_geo_grid_names)

smallRunComputing.playModel(bigRun.sub_paths,save_paths, current_run,action, current_model, intent_grid_type=1)

