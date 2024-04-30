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
# run.setActionKey("images",True)
# do_list = ["ModelC11","ModelC12","ModelC13","ModelC21","ModelC23","ModelC31","ModelC32","ModelC33"]
# do_list = ["ModelC11","ModelC12","ModelC13","ModelC21","ModelC22","ModelC23","ModelC31","ModelC32","ModelC33"]
# do_list = ["ModelC22","ModelC23","ModelC31","ModelC32","ModelC33"]
# do_list = None
do_list =None
# run.setAction(
#     {
#         "var": "nu0",
#         "start": 670e9,
#         "stop": 700e9,
#         "step": 20e9,
#         "images": True
#     }
# )
isContinuous= False

frequency_list = [86e9,230e9,345e9]
# do_list = ["ModelC22"]
# do_list = ["ModelC21","ModelC22","ModelC23"]


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
average = True

""" Intensity Grid Analysis_________________________________________________________"""
# bigRun.intensityGridAnalysis(
#     run.getAction(),
#     do_list=do_list,
#     isContinuous=isContinuous,
#     average=average
# )

""" Clean Graph Creation_________________________________________________________"""
bigRun.graphCreation(
    run.getAction(),
    do_list=do_list,
    isContinuous=False,
    average=average,
    doFullImages=False
)

""" Radial Profiles_________________________________________________________"""

# bigRun.creatRadialProfiles(
#     run.getAction(),
#     ["ModelC21","ModelC23"],
#     isContinuous=isContinuous,
#     frequency_list=frequency_list
# )

""" Blurr Intensity Grid Creation_________________________________________________________"""
#
# blurr_frequency_list = [86e9,230e9,345e9]
# blurr_kernal = [1,5,10,20]
#
# bigRun.blurrIntensityGrid(
#     run.getAction(),
#     do_list,
#     blurr_frequency_list=blurr_frequency_list,
#     blur_kernal=blurr_kernal[0]
# )
# bigRun.blurrIntensityGrid(
#     run.getAction(),
#     do_list,
#     blurr_frequency_list=blurr_frequency_list,
#     blur_kernal=blurr_kernal[1]
# )
#
# bigRun.blurrIntensityGrid(
#     run.getAction(),
#     do_list,
#     blurr_frequency_list=blurr_frequency_list,
#     blur_kernal=blurr_kernal[2]
# )
#
# bigRun.blurrIntensityGrid(
#     run.getAction(),
#     do_list,
#     blurr_frequency_list=blurr_frequency_list,
#     blur_kernal=blurr_kernal[3]
# )
# """ Blurring intensity analysis grid_________________________________________________________"""
# bigRun.blurrIntensityGridAnalysis(
#     run.getAction(),
#     do_list=do_list,
#     blurr_frequency_list=blurr_frequency_list,
#     blur_kernal=blurr_kernal[0]
# )
#
# bigRun.blurrIntensityGridAnalysis(
#     run.getAction(),
#     do_list=do_list,
#     blurr_frequency_list=blurr_frequency_list,
#     blur_kernal=blurr_kernal[1]
# )
#
# bigRun.blurrIntensityGridAnalysis(
#     run.getAction(),
#     do_list=do_list,
#     blurr_frequency_list=blurr_frequency_list,
#     blur_kernal=blurr_kernal[2]
# )
#
# bigRun.blurrIntensityGridAnalysis(
#     run.getAction(),
#     do_list=do_list,
#     blurr_frequency_list=blurr_frequency_list,
#     blur_kernal=blurr_kernal[3]
# )
#

# """ Blurr Graph Creation_________________________________________________________"""
# bigRun.blurrGraphCreation(
#     run.getAction(),
#     do_list=do_list,
#     blurr_frequency_list=blurr_frequency_list,
#     blur_kernal=blurr_kernal[0],
# )
#
# bigRun.blurrGraphCreation(
#     run.getAction(),
#     do_list=do_list,
#     blurr_frequency_list=blurr_frequency_list,
#     blur_kernal=blurr_kernal[1],
# )
#
# bigRun.blurrGraphCreation(
#     run.getAction(),
#     do_list=do_list,
#     blurr_frequency_list=blurr_frequency_list,
#     blur_kernal=blurr_kernal[2]
# )
#
# bigRun.blurrGraphCreation(
#     run.getAction(),
#     do_list=do_list,
#     blurr_frequency_list=blurr_frequency_list,
#     blur_kernal=blurr_kernal[3]
# )


