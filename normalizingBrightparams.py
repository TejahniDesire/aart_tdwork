import subprocess

from lmfit import Parameters

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
    bp["nu0"] = 230e9

    nth0_survay_points = [1.9 * 10**4, 1.3 * 10 ** 5, 8 * 10 ** 2]
    thin_total_fluxes = []
    thick_total_fluxes = []

    for i in range(len(nth0_survay_points)):
        bp["n_th0"] = nth0_survay_points[i]

        args = bigRunComputing.createIntensityArgs(bp)
        args += "--lband " + lband + " --rtray " + rtray

        subprocess.run(['python3 ' + EZPaths.aartPath + '/radialintensity.py' + args], shell=True)

        # Read created file

        fnrays = fileloading.intensityNameNoUnits(brightparams, astroModels.funckeys)
        h5f = h5py.File(fnrays, 'r')

        I0 = h5f['bghts0'][:]
        I1 = h5f['bghts1'][:]
        I2 = h5f['bghts2'][:]

        # I0_Absorb = h5f['bghts0_absorbtion'][:]
        # I1_Absorb = h5f['bghts1_absorbtion'][:]
        # I2_Absorb = h5f['bghts2_absorbtion'][:]
        Absorbtion_Image = h5f['bghts_full_absorbtion'][:]
        #
        # tau2 = h5f['tau2'][:]
        # tau1 = h5f['tau1'][:]
        # tau0 = h5f['tau0'][:]

        h5f.close()

        thin_total_fluxes += [ilp.total_jy(I0 + I1 + I2, 230e9, bp["mass"])]
        thick_total_fluxes += [ilp.total_jy(Absorbtion_Image, 230e9, bp["mass"])]





    r = np.array([2, 15, 50])
    y = full_b_func(r, mass, beta, rb_0, n_th0, p_dens).value
    params = Parameters()
    params.add('b_0', value=1)
    params.add('p_b', value=1)
    params.add('rg', value=rg_func(mass).value, vary=False)
    params.add('rb', value=rb_func(mass).value, vary=False)
    fitted_params = minimize(total_jy_normal_func, params, args=(nth0_survay_points , y), method='least_squares')
    return fitted_params.params['b_0'].value, fitted_params.params['p_b'].value


def total_jy_normal_func(params, r, y):
    rg = params['rg']
    rb = params['rb']
    b_0 = params['b_0']
    p_b = params['p_b']
    y_fit = b_0 * (r * rg / rb) ** p_b
    return y_fit - y
