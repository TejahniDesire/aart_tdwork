import matplotlib.pyplot as plt
import numpy as np
from aart_func import *
from params import *  # The file params.py contains all the relevant parameters for the simulations
from astropy import units as u
import kgeo
import image_tools as tls
import subprocess
import scipy.interpolate
from matplotlib import ticker
import importlib
from functools import partial



bp_run1 = {
    "p_mag": [-2,-1.5,-1,-0.5],
    "p_temp": [-1.2,-1,-.8,-.6],
    "p_dens": [-.7],
    "n_th0": [1.9e4],
    "t_e0": [2e11],
    "b_0": [8.131273135591028],
    "theta_b": [50.0 * (np.pi / 180)],  # NONE VARRYING________________________________________
    "mass": [(MMkg * u.kg).to(u.g).value],
    "nu0": [230e9],
    "scale_height": [.5],
    "rb_0": [2],
    "beta": [1.0],  # Legacy _____________________________________
    "r_ie": [10.0],
    "nscale": [.4]
}

bp_testRun1 = {
    "p_mag": [-2,-1.5],
    "p_temp": [-1.2,-1],
    "p_dens": [-.7],
    "n_th0": [1.9e4],
    "t_e0": [2e11],
    "b_0": [8.131273135591028],
    "theta_b": [50.0 * (np.pi / 180)],  # NONE VARRYING________________________________________
    "mass": [(MMkg * u.kg).to(u.g).value],
    "nu0": [230e9],
    "scale_height": [.5],
    "rb_0": [2],
    "beta": [1.0],  # Legacy _____________________________________
    "r_ie": [10.0],
    "nscale": [.4]
}

bp_fiducial230 = {
    "nu0": 230e9,  # 0 nu0
    "mass": (MMkg * u.kg).to(u.g).value,  # 1 mass
    "scale_height": .5,  # 2 scale_height
    "theta_b": 50.0 * (np.pi / 180),  # 3 theta_b
    "beta": 1.0,  # 4 beta
    "r_ie": 10.0,  # 5 rie
    "rb_0": 2,  # 7 rb_0
    "n_th0": 1.9e4,  # ,4e5 # 8 n_th0
    "t_e0": 2e11,  # 9 t_e0 1e12
    "p_dens": -.7,  # 10 p_dens
    "p_temp": -.84,  # 11 p_temp
    "nscale": .4  # Scale of Inoisy
}
