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


def normalize(lband,rtray,brightparams:dict):

    bp = brightparams.copy()

    fitparams = Parameters()
    fitparams.add('n_th0', value=1.3 * 10 ** 5)
    # params.add('n_th0', value=1)
    # params.add('rg', value=rg_func(mass).value, vary=False)
    fitted_params = minimize(total_jy_normal_func, fitparams, args=(lband, rtray, bp, .5), method='least_squares')
    return fitted_params.params['n_th0'].value


def total_jy_normal_func(fitparams,lband,rtray,bp,y):

    bp["n_th0"] = fitparams['n_th0'].value
    print("______________________________________ Running normalizing instance with n_th0 value of "
          + str(bp["n_th0"]) + " ______________________________________")

    args = bigRunComputing.createIntensityArgs(bp)
    args += "--lband " + lband + " --rtray " + rtray

    subprocess.run(['python3 ' + EZPaths.aartPath + '/radialintensity.py' + args], shell=True)

    # Read created file

    fnrays = fileloading.intensityNameNoUnits(bp,astroModels.funckeys)
    h5f = h5py.File(fnrays, 'r')

    I0 = h5f['bghts0'][:]
    I1 = h5f['bghts1'][:]
    I2 = h5f['bghts2'][:]

    # I0_Absorb = h5f['bghts0_absorbtion'][:]
    # I1_Absorb = h5f['bghts1_absorbtion'][:]
    # I2_Absorb = h5f['bghts2_absorbtion'][:]
    # Absorbtion_Image = h5f['bghts_full_absorbtion'][:]
    #
    # tau2 = h5f['tau2'][:]
    # tau1 = h5f['tau1'][:]
    # tau0 = h5f['tau0'][:]

    h5f.close()

    thin_total_fluxes = ilp.total_jy(I0 + I1 + I2, 230e9, bp["mass"]).value
    # thick_total_fluxes = ilp.total_jy(Absorbtion_Image, 230e9, bp["mass"])

    subprocess.run(["rm " + fnrays], shell=True)

    return y - thin_total_fluxes
