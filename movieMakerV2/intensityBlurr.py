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


def blurrIntensity(BH_mass,thin_image,absorb_image,blur_kernal):
    """

    Args:
        BH_mass: Black Hole mass in grams (without astropy units)
        thin_image: I0 + I1 + I2
        absorb_image: All passes included
        blur_kernal: In muas

    Returns:

    """
    dx = params.dx0
    one_M = ilp.rg_func(BH_mass * u.g).to(u.m)  # one Mass length unit = 1 r_g
    mass_to_uas = np.arctan(one_M.value / dBH) / muas_to_rad  # dBH is in units of meters
    # muas_blurr = 10
    rg_blurr = blur_kernal / mass_to_uas

    sig = rg_blurr / (dx * (2 * np.sqrt(
        2 * np.log(2))))  # We have 20 uas FWHM resolution. dx = uas/pixels. so 20/dx is FWHM in pixel units.
    thin_blurr_image = ndimage.gaussian_filter(thin_image, sigma=(sig, sig))
    thick_blurr_image = ndimage.gaussian_filter(absorb_image, sigma=(sig, sig))
    # print(R"$\mu a s$: " + str(mass_to_uas) + "\n"
    #       + R"$R_g$_blurr: " + str(rg_blurr) + "\n"
    #       + "sig: " + str(sig) + "\n"
    #       + "dx: " + str(dx))

    return thin_blurr_image,thick_blurr_image
