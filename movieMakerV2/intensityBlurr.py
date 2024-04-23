import subprocess

from lmfit import Parameters, minimize, fit_report
from scipy import ndimage

import EZPaths
import os

import bigRunComputing
import fileloading
from aart_func import *
import params
from params import *
from astropy import units as u
import numpy as np


def blurrIntensity(brightparams,thin_image,absorb_image):
    dx = params.dx0
    one_M = ilp.rg_func(brightparams["mass"] * u.g).to(u.m)  # one Mass length unit = 1 r_g
    mass_to_uas = np.arctan(one_M.value / dBH) / muas_to_rad  # dBH is in units of meters
    muas_blurr = 20
    # muas_blurr = 10
    rg_blurr = muas_blurr / mass_to_uas

    sig = rg_blurr / (dx * (2 * np.sqrt(
        2 * np.log(2))))  # We have 20 uas FWHM resolution. dx = uas/pixels. so 20/dx is FWHM in pixel units.
    thin_blurr_image = ndimage.gaussian_filter(thin_image, sigma=(sig, sig))
    thick_blurr_image = ndimage.gaussian_filter(absorb_image, sigma=(sig, sig))
    # print(R"$\mu a s$: " + str(mass_to_uas) + "\n"
    #       + R"$R_g$_blurr: " + str(rg_blurr) + "\n"
    #       + "sig: " + str(sig) + "\n"
    #       + "dx: " + str(dx))

    return thin_blurr_image,thick_blurr_image
