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
from final_paper_params import (
    bp_fiducial230,
    bp_fiducial345,
    bp_fiducial86,
    bp_steeperT230,
    bp_steeperT345,
    bp_steeperT86,
    bp_shallowT230,
    bp_shallowT345,
    bp_shallowT86,
    fk_fiducial
)

cmd1_args = {
    "nu0": '--nu ',
    "mass": '--mass ',
    "scale_height": '--scaleh ',
    "theta_b": '--thetab ',
    "beta": '--beta ',
    "r_ie": '--rie ',
    "rb_0": '--rb0 ',
    "n_th0": '--nth0 ',
    "t_e0": '--te0 ',
    "p_dens": '--pdens ',
    "p_temp": '--ptemp ',
    "nscale": '--nscale ',
}

cmd2_args = {
    "emodelkey": '--emodelkey ',
    "bkey": '--bkey ',
    "nnoisykey": '--nnoisykey ',
    "tnoisykey": '--tnoisykey ',
    "bnoisykey": '--bnoisykey ',
}

bps = [
    bp_fiducial230,
    bp_fiducial345,
    bp_fiducial86,
    bp_steeperT230,
    bp_steeperT345,
    bp_steeperT86,
    bp_shallowT230,
    bp_shallowT345,
    bp_shallowT86,
]
file_path = "/home/tej/Desktop/Code_Stuff/Repositories/aart_results/Final_Paper/Final_Data/"
filenames = [
    "bp_fiducial230",
    "bp_fiducial345",
    "bp_fiducial86",
    "bp_steeperT230",
    "bp_steeperT345",
    "bp_steeperT86",
    "bp_shallowT230",
    "bp_shallowT345",
    "bp_shallowT86",
    "fk_fiducial"
]
aartpath = '/home/tej/Desktop/Code_Stuff/Repositories/aart'

for i in range(len(bps)):
    print("Running: ", filenames[i])
    args = ' '

    for arg in cmd1_args:
        args = args + cmd1_args[arg] + str(bps[i][arg]) + ' '

    for arg in cmd2_args:
        args = args + cmd2_args[arg] + str(fk_fiducial[arg]) + ' '

    subprocess.run(['python3 ' + aartpath + '/radialintensity.py' + args], shell=True)

    fnrays = path + ('Intensity_a_{}_i_{}_nu_{}_mass_{}_scaleh_{}_thetab_{}_beta_{}_rie_{}_rb_{}_nth0_{}_te0_{}_pdens_{'
                     '}_ptemp_{}_nscale_{}_emkey_{}_bkey_{}_nkey_{}_tnkey_{}_bnkey_{}.h5').format(
        spin_case,
        i_case,
        "{:.5e}".format(bps[i]["nu0"]),
        "{:.5e}".format(bps[i]["mass"]),
        float(bps[i]["scale_height"]),
        "{:.3f}".format(bps[i]["theta_b"]),
        "{:.2f}".format(float(bps[i]["beta"])),
        "{:.1f}".format(float(bps[i]["r_ie"])),
        "{:.1f}".format(float(bps[i]["rb_0"])),
        "{:.1e}".format(bps[i]["n_th0"]),
        "{:.1e}".format(bps[i]["t_e0"]),
        float(bps[i]["p_dens"]),
        float(bps[i]["p_temp"]),
        "{:.1f}".format(bps[i]["nscale"]),
        fk_fiducial["emodelkey"],
        fk_fiducial["bkey"],
        fk_fiducial["nnoisykey"],
        fk_fiducial["tnoisykey"],
        fk_fiducial["bnoisykey"]
    )

    print("Reading file: ", fnrays)

    h5f = h5py.File(fnrays, 'r')

    I0 = h5f['bghts0'][:]  # This implies I0 is 1 pass
    I1 = h5f['bghts1'][:]
    I2 = h5f['bghts2'][:]

    I2_Absorb = h5f['bghts2_absorbtion'][:]
    I1_Absorb = h5f['bghts1_absorbtion'][:]
    I0_Absorb = h5f['bghts0_absorbtion'][:]
    Absorbtion_Image = h5f['bghts_full_absorbtion'][:]
    full_profiles0 = h5f['full_profiles0'][:]
    full_profiles1 = h5f['full_profiles1'][:]
    full_profiles2 = h5f['full_profiles2'][:]
    full_profiles_unit = h5f['full_profiles_unit'][:]

    h5f.close()

    h5f = h5py.File(file_path + filenames[i], 'w')

    # h5f.create_dataset('I0', data=I0)
    # h5f.create_dataset('I1', data=I1)
    # h5f.create_dataset('I2', data=I2)
    #
    # h5f.create_dataset('I0_Absorb', data=I0_Absorb)
    # h5f.create_dataset('I1_Absorb', data=I1_Absorb)
    # h5f.create_dataset('I2_Absorb', data=I2_Absorb)
    # h5f.create_dataset('Absorbtion_Image', data=Absorbtion_Image)

    h5f.create_dataset('full_profiles0', data=full_profiles0)
    h5f.create_dataset('full_profiles1', data=full_profiles1)
    h5f.create_dataset('full_profiles2', data=full_profiles2)

    h5f.close()

    print("File {} created".format(file_path + filenames[i]))








