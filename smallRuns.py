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


run = runDataClass.soloRun2
do_list = None
isContinuous= False


run.setisNormalized(True)
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
#     isContinuous=isContinuous,
#     frequency_list=None
# )

""" Intensity Grid Analysis_________________________________________________________"""
bigRun.intensityGridAnalysis(
    run.getAction(),
    do_list=do_list,
    isContinuous=isContinuous,
    average=False
)
""" Clean Graph Creation_________________________________________________________"""
# bigRun.graphCreation(
#     run.getAction(),
#     do_list=do_list,
#     isContinuous=False
# )
