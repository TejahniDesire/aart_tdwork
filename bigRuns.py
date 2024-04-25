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
do_list = ["ModelC11","ModelC12","ModelC13","ModelC21","ModelC22","ModelC23","ModelC31","ModelC32","ModelC33"]
# do_list = ["ModelC22","ModelC23","ModelC31","ModelC32","ModelC33"]
# do_list = None

bigRun = classRunComputing.BigRuns(
    run.getRunName(),
    run.getBrightparams(),
    run.getBPVarNames(),
    run.getGeoGrid(),
    run.getGeoGridNames(),
    normalized_brightparams=run.getIsNormalized(),
)
""" Geo Model_________________________________________________________"""

# bigRun.createGeoGrid()

""" Intensity Grid Creation_________________________________________________________"""
#
# bigRun.creatIntensityGrid(
#     run.getAction(),
#     do_list=do_list,
#     isContinuous=True
# )

# bigRun.blurrIntensityGrid(
#     run.getAction(),
#     do_list
# )

""" Intensity Grid Analysis_________________________________________________________"""
bigRun.intensityGridAnalysis(
    run.getAction(),
    do_list=do_list,
    isContinuous=True
)
""" Clean Graph Creation_________________________________________________________"""
bigRun.graphCreation(
    run.getAction(),
    do_list=do_list,
    isContinuous=True
)

""" Blurring intensity grid_________________________________________________________"""
# bigRun.blurrIntensityGridAnalysis(
#     run.getAction(),
#     do_list=["ModelA11","ModelB31"],
#     isContinuous=True
# )

""" Blurr Graph Creation_________________________________________________________"""
#
# bigRun.blurrGraphCreation(
#     run.getAction()
# )
""" _________________________________________________________"""

