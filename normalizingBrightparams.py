import subprocess
import EZPaths
import os
from aart_func import *
import params
from params import *
import importlib
import astroModels


def normalize(brightparams:dict):
    bp = brightparams.copy()
    bp["nu0"] = 230e9

    nth0_survay_points = [1.9 * 10**4, 1.3 * 10 ** 5, 8 * 10 ** 2]
    # y =
    #
    # r = np.array([2, 15, 50])
    # y = full_b_func(r, mass, beta, rb_0, n_th0, p_dens).value
    # params = Parameters()
    # params.add('b_0', value=1)
    # params.add('p_b', value=1)
    # params.add('rg', value=rg_func(mass).value, vary=False)
    # params.add('rb', value=rb_func(mass).value, vary=False)
    # fitted_params = minimize(b_simple_func, params, args=(r, y), method='least_squares')
    # return fitted_params.params['b_0'].value, fitted_params.params['p_b'].value
