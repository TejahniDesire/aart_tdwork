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


run = runDataClass.run2
do_list = ["ModelB21","ModelB22","ModelB23"]

bigRun = classRunComputing.BigRuns(
    run.getRunName(),
    run.getBrightparams(),
    run.getBPVarNames(),
    run.getGeoGrid(),
    run.getGeoGridNames(),
    normalized_brightparams=run.getIsNormalized(),
    isContinuous=run.getIsContinuous()
)
""" Geo Model_________________________________________________________"""

# bigRun.createGeoGrid()

""" Intensity Grid Creation_________________________________________________________"""

# bigRun.creatIntensityGrid(
#     run.getAction(),
#     do_list=do_list
# )

# bigRun.blurrIntensityGrid(
#     run.getAction()
# )

""" Intensity Grid Analysis_________________________________________________________"""
# bigRun.intensityGridAnalysis(
#     run.getAction(),
#     do_list=do_list
# )
""" Clean Graph Creation_________________________________________________________"""
bigRun.graphCreation(
    run.getAction()
)

""" Blurring intensity grid_________________________________________________________"""
# bigRun.blurrIntensityGridAnalysis(
#     run.getAction()
# )

""" Blurr Graph Creation_________________________________________________________"""

# bigRun.blurrGraphCreation(
#     run.getAction()
# )
""" _________________________________________________________"""

