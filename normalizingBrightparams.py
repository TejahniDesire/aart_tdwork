import subprocess

from lmfit import Parameters, minimize, fit_report

import EZPaths
import os

import bigRunComputing
import fileloading
from aart_func import *
import params
from params import *
import importlib
import astroModels
from movieMakerV2 import  intensityBlurr


def normalize(lband,rtray,brightparams:dict):

    bp = brightparams.copy()
    bp["nu0"] = 230e9
    fitparams = Parameters()
    fitparams.add('n_th0', value=1.3 * 10 ** 5,min=0)
    fitted_params = minimize(total_jy_normal_func, fitparams, args=(lband, rtray, bp, .5), method='least_squares')
    return fitted_params.params['n_th0'].value


def total_jy_normal_func(fitparams,lband,rtray,bp,y):

    bp["n_th0"] = fitparams['n_th0'].value
    print("______________________________________ Running normalizing instance with n_th0 value of "
          + str(bp["n_th0"]) + " ______________________________________")

    thin_total_flux,thick_total_flux = totalIntensity230Point(lband,rtray,bp,True)
    print("Total Thin Model Flux of ",thin_total_flux)
    print("Total Full Model Flux of ",thick_total_flux)

    return y - thick_total_flux


def totalIntensity230Point(lband,rtray,brightparams:dict,already230=False,blurr_policy=False,blur_kernal=None):
    if already230:
        bp = brightparams
    else:
        bp = brightparams.copy()
        bp["nu0"] = 230e9

    args = bigRunComputing.createIntensityArgs(bp)
    args += "--lband " + lband + " --rtray " + rtray

    subprocess.run(['python3 ' + EZPaths.aartPath + '/radialintensity.py' + args], shell=True)

    # Read created file

    fnrays = fileloading.intensityNameNoUnits(bp, astroModels.funckeys)
    h5f = h5py.File(fnrays, 'r')

    I0 = h5f['bghts0'][:]
    I1 = h5f['bghts1'][:]
    I2 = h5f['bghts2'][:]

    # I0_Absorb = h5f['bghts0_absorbtion'][:]
    # I1_Absorb = h5f['bghts1_absorbtion'][:]
    # I2_Absorb = h5f['bghts2_absorbtion'][:]
    absorb_image = h5f['bghts_full_absorbtion'][:]
    thin_image = I0 + I1 + I2
    #
    # tau2 = h5f['tau2'][:]
    # tau1 = h5f['tau1'][:]
    # tau0 = h5f['tau0'][:]

    if blurr_policy:
        thin_image, absorb_image = intensityBlurr.blurrIntensity(brightparams["mass"],thin_image,absorb_image,blur_kernal)

    h5f.close()

    subprocess.run(["rm " + fnrays], shell=True)

    thin_total_flux = ilp.total_jy(thin_image,230e9,bp["mass"]).value
    thick_total_flux = ilp.total_jy(absorb_image,230e9,bp["mass"]).value

    return thin_total_flux,thick_total_flux
