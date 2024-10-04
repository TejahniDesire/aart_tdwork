import sys
import subprocess
import EZPaths
import os
import importlib
from functools import partial 
from importlib import reload

aartpath = '/home/td6241/repositories/aart' #insert path to aart repo
sys.path.append(aartpath)

def loadGeoModel(current_model:str, run:str):
    print("\nLoading " + current_model + "\n")
    paramFile = EZPaths.aartPath + "/params.py"
    model_is_loaded = os.path.exists(paramFile)
    if model_is_loaded:
        subprocess.run(["rm " + paramFile], shell=True)

    loading_model = EZPaths.modelRunsDir + run + "/geoModels/" + current_model + ".py"
    loaded_model = EZPaths.aartPath + "/" + "params.py"
    cmd = "cp " + loading_model + " " + loaded_model
    subprocess.run([cmd], shell=True)
    print("file {} Loaded as {}".format(loading_model,loaded_model))
    print("Remember to call 'importlib.reload(params)' to reload parameters")
    return partial(exec,"import params; importlib.reload( params ); from params import *")

def reloadVariables():
    return "import params; importlib.reload( params ); from params import *"