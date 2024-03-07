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


fk_fiducial = {
    "emodelkey": 0,  # emodelkey Emission Model choice, 0 = thermal ultrarelativistic, 1 = power law
    "bkey": 2,  # bkey
    "nnoisykey": 0,  # nnoisykey Inoisy density. 0 = no noise, 1 = noise
    "tnoisykey": 0,  # tnoisykey Inoisy temperature
    "bnoisykey": 0  # bnoisykey Inoisy magnetic field
}
# kw_bv_0 = 8.131273135591028
# kw_p_bv = -.3

"""None Units_________________________________________________________________________________________________"""

# 230e9________________________________________________________________________________________________
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

# Model 2
bp_steeperT230 = {
    "nu0": 230e9,  # 0 nu0
    "mass": (MMkg * u.kg).to(u.g).value,  # 1 mass
    "scale_height": .5,  # 2 scale_height
    "theta_b": 50.0 * (np.pi / 180),  # 3 theta_b
    "beta": 1.0,  # 4 beta
    "r_ie": 10.0,  # 5 rie
    "rb_0": 2,  # 7 rb_0
    "n_th0": 1.3e5,  # ,4e5 # 8 n_th0
    "t_e0": 2e11,  # 9 t_e0 1e12
    "p_dens": -.7,  # 10 p_dens
    "p_temp": -1.6,  # 11 p_temp
    "nscale": .4  # Scale of Inoisy
}

# Model 3
bp_shallowT230 = {
    "nu0": 230e9,  # 0 nu0
    "mass": (MMkg * u.kg).to(u.g).value,  # 1 mass
    "scale_height": .5,  # 2 scale_height
    "theta_b": 50.0 * (np.pi / 180),  # 3 theta_b
    "beta": 1.0,  # 4 beta
    "r_ie": 10.0,  # 5 rie
    "rb_0": 2,  # 7 rb_0
    "n_th0": 8e2,  # ,4e5 # 8 n_th0
    "t_e0": 2e11,  # 9 t_e0 1e12
    "p_dens": -.7,  # 10 p_dens
    "p_temp": -.3,  # 11 p_temp
    "nscale": .4  # Scale of Inoisy
}


# 345e9________________________________________________________________________________________________
bp_fiducial345 = {
    "nu0": 345e9,  # 0 nu0
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

# Model 2
bp_steeperT345 = {
    "nu0": 345e9,  # 0 nu0
    "mass": (MMkg * u.kg).to(u.g).value,  # 1 mass
    "scale_height": .5,  # 2 scale_height
    "theta_b": 50.0 * (np.pi / 180),  # 3 theta_b
    "beta": 1.0,  # 4 beta
    "r_ie": 10.0,  # 5 rie
    "rb_0": 2,  # 7 rb_0
    "n_th0": 1.3e5,  # ,4e5 # 8 n_th0
    "t_e0": 2e11,  # 9 t_e0 1e12
    "p_dens": -.7,  # 10 p_dens
    "p_temp": -1.6,  # 11 p_temp
    "nscale": .4  # Scale of Inoisy
}

# Model 3
bp_shallowT345 = {
    "nu0": 345e9,  # 0 nu0
    "mass": (MMkg * u.kg).to(u.g).value,  # 1 mass
    "scale_height": .5,  # 2 scale_height
    "theta_b": 50.0 * (np.pi / 180),  # 3 theta_b
    "beta": 1.0,  # 4 beta
    "r_ie": 10.0,  # 5 rie
    "rb_0": 2,  # 7 rb_0
    "n_th0": 8e2,  # ,4e5 # 8 n_th0
    "t_e0": 2e11,  # 9 t_e0 1e12
    "p_dens": -.7,  # 10 p_dens
    "p_temp": -.3,  # 11 p_temp
    "nscale": .4  # Scale of Inoisy
}


# 86e9________________________________________________________________________________________________
bp_fiducial86 = {
    "nu0": 86e9,  # 0 nu0
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

# Model 2
bp_steeperT86 = {
    "nu0": 86e9,  # 0 nu0
    "mass": (MMkg * u.kg).to(u.g).value,  # 1 mass
    "scale_height": .5,  # 2 scale_height
    "theta_b": 50.0 * (np.pi / 180),  # 3 theta_b
    "beta": 1.0,  # 4 beta
    "r_ie": 10.0,  # 5 rie
    "rb_0": 2,  # 7 rb_0
    "n_th0": 1.3e5,  # ,4e5 # 8 n_th0
    "t_e0": 2e11,  # 9 t_e0 1e12
    "p_dens": -.7,  # 10 p_dens
    "p_temp": -1.6,  # 11 p_temp
    "nscale": .4  # Scale of Inoisy
}

# Model 3
bp_shallowT86 = {
    "nu0": 86e9,  # 0 nu0
    "mass": (MMkg * u.kg).to(u.g).value,  # 1 mass
    "scale_height": .5,  # 2 scale_height
    "theta_b": 50.0 * (np.pi / 180),  # 3 theta_b
    "beta": 1.0,  # 4 beta
    "r_ie": 10.0,  # 5 rie
    "rb_0": 2,  # 7 rb_0
    "n_th0": 8e2,  # ,4e5 # 8 n_th0
    "t_e0": 2e11,  # 9 t_e0 1e12
    "p_dens": -.7,  # 10 p_dens
    "p_temp": -.3,  # 11 p_temp
    "nscale": .4  # Scale of Inoisy
}
"""Units___________________________________________________________________________________________"""

bp_fiducial_unit = {
    "nu0": 230e9*ilp.Hz,  # 0 nu0
    "mass": (MMkg * u.kg).to(u.g),  # 1 mass
    "scale_height": .5,  # 2 scale_height
    "theta_b": 50.0 * (np.pi / 180)*ilp.rads,  # 3 theta_b
    "beta": 1.0,  # 4 beta
    "r_ie": 10.0,  # 5 rie
    "rb_0": 2,  # 7 rb_0
    "n_th0": 1.9e4*ilp.cmcubed,  # ,4e5 # 8 n_th0
    "t_e0": 2e11*ilp.kelv,  # 9 t_e0 1e12
    "p_dens": -.7,  # 10 p_dens
    "p_temp": -.84,  # 11 p_temp
    "nscale": .4  # Scale of Inoisy
}
# kw_bv_0 = 8.131273135591028
# kw_p_bv = -.3

# Model 2
bp_steeperT_unit = {
    "nu0": 230e9*ilp.Hz,  # 0 nu0
    "mass": (MMkg * u.kg).to(u.g),  # 1 mass
    "scale_height": .5,  # 2 scale_height
    "theta_b": 50.0 * (np.pi / 180)*ilp.rads,  # 3 theta_b
    "beta": 1.0,  # 4 beta
    "r_ie": 10.0,  # 5 rie
    "rb_0": 2,  # 7 rb_0
    "n_th0": 1.3e5*ilp.cmcubed,  # ,4e5 # 8 n_th0
    "t_e0": 2e11*ilp.kelv,  # 9 t_e0 1e12
    "p_dens": -.7,  # 10 p_dens
    "p_temp": -1.6,  # 11 p_temp
    "nscale": .4  # Scale of Inoisy
}

# Model 3
bp_shallowT_unit = {
    "nu0": 230e9*ilp.Hz,  # 0 nu0
    "mass": (MMkg * u.kg).to(u.g),  # 1 mass
    "scale_height": .5,  # 2 scale_height
    "theta_b": 50.0 * (np.pi / 180)*ilp.rads,  # 3 theta_b
    "beta": 1.0,  # 4 beta
    "r_ie": 10.0,  # 5 rie
    "rb_0": 2,  # 7 rb_0
    "n_th0": 8e2*ilp.cmcubed,  # ,4e5 # 8 n_th0
    "t_e0": 2e11*ilp.kelv,  # 9 t_e0 1e12
    "p_dens": -.7,  # 10 p_dens
    "p_temp": -.3,  # 11 p_temp
    "nscale": .4  # Scale of Inoisy
}

# 345e9____________________________________________________________________________________________________

bp_fiducial_345unit = {
    "nu0": 345e9*ilp.Hz,  # 0 nu0
    "mass": (MMkg * u.kg).to(u.g),  # 1 mass
    "scale_height": .5,  # 2 scale_height
    "theta_b": 50.0 * (np.pi / 180)*ilp.rads,  # 3 theta_b
    "beta": 1.0,  # 4 beta
    "r_ie": 10.0,  # 5 rie
    "rb_0": 2,  # 7 rb_0
    "n_th0": 1.9e4*ilp.cmcubed,  # ,4e5 # 8 n_th0
    "t_e0": 2e11*ilp.kelv,  # 9 t_e0 1e12
    "p_dens": -.7,  # 10 p_dens
    "p_temp": -.84,  # 11 p_temp
    "nscale": .4  # Scale of Inoisy
}

# Model 2
bp_steeperT_345unit = {
    "nu0": 345e9*ilp.Hz,  # 0 nu0
    "mass": (MMkg * u.kg).to(u.g),  # 1 mass
    "scale_height": .5,  # 2 scale_height
    "theta_b": 50.0 * (np.pi / 180)*ilp.rads,  # 3 theta_b
    "beta": 1.0,  # 4 beta
    "r_ie": 10.0,  # 5 rie
    "rb_0": 2,  # 7 rb_0
    "n_th0": 1.3e5*ilp.cmcubed,  # ,4e5 # 8 n_th0
    "t_e0": 2e11*ilp.kelv,  # 9 t_e0 1e12
    "p_dens": -.7,  # 10 p_dens
    "p_temp": -1.6,  # 11 p_temp
    "nscale": .4  # Scale of Inoisy
}

# Model 3
bp_shallowT_345unit = {
    "nu0": 345e9*ilp.Hz,  # 0 nu0
    "mass": (MMkg * u.kg).to(u.g),  # 1 mass
    "scale_height": .5,  # 2 scale_height
    "theta_b": 50.0 * (np.pi / 180)*ilp.rads,  # 3 theta_b
    "beta": 1.0,  # 4 beta
    "r_ie": 10.0,  # 5 rie
    "rb_0": 2,  # 7 rb_0
    "n_th0": 8e2*ilp.cmcubed,  # ,4e5 # 8 n_th0
    "t_e0": 2e11*ilp.kelv,  # 9 t_e0 1e12
    "p_dens": -.7,  # 10 p_dens
    "p_temp": -.3,  # 11 p_temp
    "nscale": .4  # Scale of Inoisy
}

# 86e9____________________________________________________________________________________________________

bp_fiducial_86unit = {
    "nu0": 86e9*ilp.Hz,  # 0 nu0
    "mass": (MMkg * u.kg).to(u.g),  # 1 mass
    "scale_height": .5,  # 2 scale_height
    "theta_b": 50.0 * (np.pi / 180)*ilp.rads,  # 3 theta_b
    "beta": 1.0,  # 4 beta
    "r_ie": 10.0,  # 5 rie
    "rb_0": 2,  # 7 rb_0
    "n_th0": 1.9e4*ilp.cmcubed,  # ,4e5 # 8 n_th0
    "t_e0": 2e11*ilp.kelv,  # 9 t_e0 1e12
    "p_dens": -.7,  # 10 p_dens
    "p_temp": -.84,  # 11 p_temp
    "nscale": .4  # Scale of Inoisy
}
# kw_bv_0 = 8.131273135591028
# kw_p_bv = -.3

# Model 2
bp_steeperT_86unit = {
    "nu0": 86e9*ilp.Hz,  # 0 nu0
    "mass": (MMkg * u.kg).to(u.g),  # 1 mass
    "scale_height": .5,  # 2 scale_height
    "theta_b": 50.0 * (np.pi / 180)*ilp.rads,  # 3 theta_b
    "beta": 1.0,  # 4 beta
    "r_ie": 10.0,  # 5 rie
    "rb_0": 2,  # 7 rb_0
    "n_th0": 1.3e5*ilp.cmcubed,  # ,4e5 # 8 n_th0
    "t_e0": 2e11*ilp.kelv,  # 9 t_e0 1e12
    "p_dens": -.7,  # 10 p_dens
    "p_temp": -1.6,  # 11 p_temp
    "nscale": .4  # Scale of Inoisy
}

# Model 3
bp_shallowT_86unit = {
    "nu0": 86e9*ilp.Hz,  # 0 nu0
    "mass": (MMkg * u.kg).to(u.g),  # 1 mass
    "scale_height": .5,  # 2 scale_height
    "theta_b": 50.0 * (np.pi / 180)*ilp.rads,  # 3 theta_b
    "beta": 1.0,  # 4 beta
    "r_ie": 10.0,  # 5 rie
    "rb_0": 2,  # 7 rb_0
    "n_th0": 8e2*ilp.cmcubed,  # ,4e5 # 8 n_th0
    "t_e0": 2e11*ilp.kelv,  # 9 t_e0 1e12
    "p_dens": -.7,  # 10 p_dens
    "p_temp": -.3,  # 11 p_temp
    "nscale": .4  # Scale of Inoisy
}



