import sys
import subprocess

aartpath = '/home/td6241/repositories/aart'  # insert path to aart repo
sys.path.append(aartpath)

import fileloading
import astroModels
import image_tools

from matplotlib import ticker
from aart_func import *
from params import *  # The file params.py contains all the relevant parameters for the simulations
from astropy import units as u
import tools as tls
import numpy as np

import subprocess
from aart_func import *
import params
from params import *  # The file params.py contains all the relevant parameters for the simulations
import image_tools as tls
import numpy as np
from movieMakerV2 import intensityBlurr

# TODO: Add counter to current value in image
# TODO: Make sig fig of counter vary with cmd_arg
# TODO: Add bar every certain amount of images
# sudo apt install ffmpeg


'''CMD ARGS'''
# Which var, Start, Stop, amount of images
parser = argparse.ArgumentParser(description='Movies', formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('var', help=
''' 
0: nu0 (Hertz)
1: mass (grams)
2: scale height (rg)
3: theta b (radians)
4: beta (dimensionless)
5: Rie (dimensionless)
6: Bchoi [DO NOT VARY]
7: rb (rg)
8: nth0 (1/cm^3)
9: te0 (Kelvin)
10: pdens (dimensionless)
11: ptemp (dimensionless)
''',
                    type=int)
parser.add_argument('start', type=float)
parser.add_argument('stop', help='Inclusive Stop', type=float)
parser.add_argument("step_size", type=float)
parser.add_argument("power_val", type=float,
                    help="Colorbar power norm value. Lower for anticipated higher brightness range")
args = parser.parse_args()
action = [
    (args.var),
    (args.start),
    (args.stop),
    (args.step_size),
    (args.power_val)
]
action = {
    'var': args.var,
    'start': args.start,
    'stop': args.stop,
    'step': args.step_size,
    'power': args.power_val
}


def curve_params(varphi, rho):
    """calculate Appendix B parameters for a curve rho(varphi)
       assume varphis are equally spaced!!!"""

    # spacing in varphi  
    dvarphi = varphi[-1] - varphi[-2]

    # area
    area = np.trapz(0.5 * rho ** 2, dx=dvarphi)

    # centroid
    mux = np.trapz((rho ** 3 * np.cos(varphi)) / (3 * area), dx=dvarphi)
    muy = np.trapz((rho ** 3 * np.sin(varphi)) / (3 * area), dx=dvarphi)

    # second moment
    Sxx = np.trapz((rho ** 4 * np.cos(varphi) ** 2) / (4 * area), dx=dvarphi) - mux ** 2
    Syy = np.trapz((rho ** 4 * np.sin(varphi) ** 2) / (4 * area), dx=dvarphi) - muy ** 2
    Sxy = np.trapz((rho ** 4 * np.sin(varphi) * np.cos(varphi)) / (4 * area), dx=dvarphi) - mux * muy

    # diagonalize 2nd moment matrix
    D = np.sqrt((Sxx - Syy) ** 2 + 4 * Sxy * Sxy)
    a = np.sqrt(2 * (Sxx + Syy + D))
    b = np.sqrt(2 * (Sxx + Syy - D))

    # radius, eccentricity, position angle
    r = np.sqrt(0.5 * (a ** 2 + b ** 2))
    e = np.sqrt(1 - b ** 2 / a ** 2)
    chi = 0.5 * np.arcsin(2 * Sxy / D)

    return (area, mux, muy, r, e, chi)


'''Reading of the lensing bands----------------------------------'''
lband = "/scratch/gpfs/td6241/aart/rawResults/LensingBands_a_%s_i_%s.h5" % (spin_case, i_case)

print("Reading file: ", lband)

h5f = h5py.File(lband, 'r')

# Points for the boundary of the BH shadow
alpha_critc = h5f['alpha'][:]
beta_critc = h5f['beta'][:]

# The concave hulls for the lensing bands
hull_0i = h5f['hull_0i'][:]
hull_0e = h5f['hull_0e'][:]
hull_1i = h5f['hull_1i'][:]
hull_1e = h5f['hull_1e'][:]
hull_2i = h5f['hull_2i'][:]
hull_2e = h5f['hull_2e'][:]

# The grid points for each lensing band
supergrid0 = h5f['grid0'][:]
N0 = int(h5f["N0"][0])
mask0 = h5f['mask0'][:]
lim0 = int(h5f["lim0"][0])
supergrid1 = h5f['grid1'][:]
N1 = int(h5f["N1"][0])
mask1 = h5f['mask1'][:]
lim1 = int(h5f["lim1"][0])
supergrid2 = h5f['grid2'][:]
N2 = int(h5f["N2"][0])
mask2 = h5f['mask2'][:]
lim2 = int(h5f["lim2"][0])

h5f.close()

'''Reading Analytical Ray-tracing----------------------------------'''
rtray = "/scratch/gpfs/td6241/aart/rawResults/Rays_a_%s_i_%s.h5" % (spin_case, i_case)

print("Reading file: ", rtray)

h5f = h5py.File(rtray, 'r')

rs0 = h5f['rs0'][:]
sign0 = h5f['sign0'][:]
t0 = h5f['t0'][:]
phi0 = h5f['phi0'][:]

rs1 = h5f['rs1'][:]
sign1 = h5f['sign1'][:]
t1 = h5f['t1'][:]
phi1 = h5f['phi1'][:]

rs2 = h5f['rs2'][:]
sign2 = h5f['sign2'][:]
t2 = h5f['t2'][:]
phi2 = h5f['phi2'][:]

h5f.close()

'''Computing images----------------------------------'''
args = ' '
cmd_args = [
    'nu',
    'mass',
    'scaleh',
    'thetab',
    'beta',
    'Rie',
    'Bchoi',
    'rb',
    'nth0',
    'te0',
    'pdens',
    'ptemp'
]
# brightparams = [
# 	230e9, # nu0
# 	(MMkg * u.kg).to(u.g).value, # mass
# 	.5, # scale_height
# 	50.0 * (np.pi / 180), # theta_b
# 	1.0, # beta
# 	10.0, # Rie
# 	0, # Bchoice
# 	50.0, # rb
# 	1.23e4, # n_th0
# 	8.1e9, # t_e0
# 	-.7, # p_dens
# 	-.84 # p_temp
# ]

# #  LOWER PTEMP
# brightparams = [
# 	230e9, # nu0
# 	(MMkg * u.kg).to(u.g).value, # mass
# 	.5, # scale_height
# 	50.0 * (np.pi / 180), # theta_b
# 	1.0, # beta
# 	10.0, # Rie
# 	0, # Bchoice
# 	2, # rb
# 	2.9726e+05, # n_th0
# 	1.2428e+11, # t_e0
# 	-.7, # p_dens
# 	-1.6 # p_temp
# ]

# # HIGHER PTEMP
brightparamsLEG = [
    230e9,  # nu0
    (MMkg * u.kg).to(u.g).value,  # mass
    .5,  # scale_height
    50.0 * (np.pi / 180),  # theta_b
    1.0,  # beta
    10.0,  # Rie
    0,  # Bchoice
    2,  # rb
    8248.16518597441,  # n_th0
    1.2428e+11,  # t_e0
    -.7,  # p_dens
    -1  # p_temp
]
brightparams = {
    "p_mag": -1.5,
    "p_temp": -1,
    "p_dens": -.7,
    "n_th0": 1.9e4,
    "t_e0": 1e11,
    "b_0": 5,
    "theta_b": 50.0 * (np.pi / 180),  # NONE VARRYING________________________________________
    "mass": (MMkg * u.kg).to(u.g).value,
    "nu0": 230e9,
    "scale_height": .5,
    "rb_0": 5,
    "beta": 1.0,  # Legacy _____________________________________
    "r_ie": 10.0,
    "nscale": .4
}

# [Var, Unit Magnitude, Units]
label = np.zeros([len(brightparamsLEG), 3], dtype=object)
label[:, 0] = [
    r"$\nu= $",
    'BlkHole Mass= ',
    'Scale Height= ',
    r'$\theta_b= $',
    r'$\beta= $',
    r'$R_ie= $',
    r'BChoi= $',
    r'$R_b= $',
    r'$n_{th,0}= $',
    r'$T_{e,0}= $',
    r'$p_{dens}= $',
    r'$p_{temp}= $'
]
label[:, 1] = [
    1e9,  # GHz
    1e9 * (1.989e33),  # Billion Solar Masses
    1,  # Rg
    1,  # Rads
    1,
    1,
    1,
    1,  # Rg
    1 / 1e6,  # 1/m^3
    1e9,  # GK
    1,
    1
]
label[:, 2] = [
    'GHz',  # GHz
    r'Billion $M_{\odot}$',  # Billion Solar Masses
    r'$R_g$',  # Rg
    'Rads',  # Rads
    '',
    '',
    '',
    r'$R_g$',  # Rg
    r'$1/m^{3}$',  # 1/m^3
    'GK',  # GK
    '',
    ''
]

num_iterations = (action[2] - action[1]) / action[3]
doth5_files = []  # .h5 files
names_to_delete2 = []  # images
b = 0
thetapointsamount = 100

x_variable = []

jansky_variable_total = []
jansky_variable_n0 = []
jansky_variable_n1 = []
jansky_variable_n2 = []

ring_radii_n0 = []
ring_radii_n1 = []
ring_radii_n2 = []

# Ignore array[0,:]
ring_radii_n0_array = np.zeros(thetapointsamount)
ring_radii_n1_array = np.zeros(thetapointsamount)
ring_radii_n2_array = np.zeros(thetapointsamount)

num_iterations = int((action["stop"] - action["start"]) / action["step"])
num_of_theta_points = image_tools.num_of_theta_points
"""GRAPHS________________________________________________________________________________________________________"""

x_variable = np.zeros(num_iterations)  # counter for independant variable

# flux_________________________________
janksys_thick = np.ndarray([num_iterations, 4])  # [I0, I1, I2, FullImage]
janksys_thin = np.ndarray([num_iterations, 4])  # [I0, I1, I2, total]

# Average Radius over all Theta_________________________________
mean_radii_Thin = np.ndarray([num_iterations, 4])  # [I0, I1, I2, FullImage]
mean_radii_Thick = np.ndarray([num_iterations, 4])  # [I0, I1, I2, FullImage]

# First layer of Radius as a function of theta_________________________________
radii_Full_Thin = np.zeros(num_of_theta_points)
radii_I0_Thin = np.zeros(num_of_theta_points)
radii_I1_Thin = np.zeros(num_of_theta_points)
radii_I2_Thin = np.zeros(num_of_theta_points)

radii_FullAbsorption_Thick = np.zeros(num_of_theta_points)
radii_I0_Thick = np.zeros(num_of_theta_points)
radii_I1_Thick = np.zeros(num_of_theta_points)
radii_I2_Thick = np.zeros(num_of_theta_points)

# Optical Depth_________________________________
mean_optical_depth_I0 = np.zeros(num_iterations)
mean_optical_depth_I1 = np.zeros(num_iterations)
mean_optical_depth_I2 = np.zeros(num_iterations)

for i in range(num_iterations):
    brightparams[action["var"]] = action["start"] + i * action["step"]
    current_freqeuncy = brightparams[action["var"]]

    x_variable[i] = current_freqeuncy

    rtray = fileloading.intensityNameNoUnits(brightparams, astroModels.funckeys)

    doth5_files += [rtray]

    print("Reading file: ", rtray)

    h5f = h5py.File(rtray, 'r')

    I0 = h5f['bghts0'][:]
    I1 = h5f['bghts1'][:]
    I2 = h5f['bghts2'][:]

    I0_Absorb = h5f['bghts0_absorbtion'][:]
    I1_Absorb = h5f['bghts1_absorbtion'][:]
    I2_Absorb = h5f['bghts2_absorbtion'][:]

    Absorbtion_Image = h5f['bghts_full_absorbtion'][:]

    h5f.close()

    # Thin Radii Calcs----------------------------------------------------------------------------------------------

    # radii_I0_Thin_i, theta = tls.radii_of_thetaV2(I0, params.dx0, average=average)
    # radii_I1_Thin_i, theta = tls.radii_of_thetaV2(I1, params.dx0, average=average)
    # radii_I2_Thin_i, theta = tls.radii_of_thetaV2(I2, params.dx0, average=average)
    # radii_Full_Thick_i, theta = tls.radii_of_thetaV2(I0 + I1 + I2, params.dx0, average=average)
    #
    # # profs = scipy.ndimage.convolve1d(profs, np.ones(navg_ang), axis=0) / navg_ang
    #
    # r0_thin = tls.curve_params(theta, radii_I0_Thin_i)
    # r1_thin = tls.curve_params(theta, radii_I1_Thin_i)
    # r2_thin = tls.curve_params(theta, radii_I2_Thin_i)
    # full_thin = tls.curve_params(theta, radii_Full_Thick_i)
    #
    # mean_radii_Thin[i, 0] = r0_thin
    # mean_radii_Thin[i, 1] = r1_thin
    # mean_radii_Thin[i, 2] = r2_thin
    # mean_radii_Thin[i, 3] = full_thin
    #
    # # radii_I0_Thin = np.vstack((radii_I0_Thin, radii_I0_Thin_i))
    # # radii_I1_Thin = np.vstack((radii_I1_Thin, radii_I1_Thin_i))
    # # radii_I2_Thin = np.vstack((radii_I2_Thin, radii_I2_Thin_i))
    # # radii_Full_Thin = np.vstack((radii_Full_Thin, radii_Full_Thick_i))

    # Thick Radii Calcs---------------------------------------------------------------------------------------------
    radii_I0_Thick_i, theta = tls.radii_of_thetaV2(I0_Absorb, params.dx0)
    radii_I1_Thick_i, theta = tls.radii_of_thetaV2(I1_Absorb, params.dx0)
    radii_I2_Thick_i, theta = tls.radii_of_thetaV2(I2_Absorb, params.dx0)
    radii_FullAbsorption_Thick_i, theta = tls.radii_of_thetaV2(Absorbtion_Image, params.dx0)

    r0_thick = tls.curve_params(theta, radii_I0_Thick_i)
    r1_thick = tls.curve_params(theta, radii_I1_Thick_i)
    r2_thick = tls.curve_params(theta, radii_I2_Thick_i)
    full_thick = tls.curve_params(theta, radii_FullAbsorption_Thick_i)

    mean_radii_Thick[i, 0] = r0_thick
    mean_radii_Thick[i, 1] = r1_thick
    mean_radii_Thick[i, 2] = r2_thick
    mean_radii_Thick[i, 3] = full_thick

    radii_I0_Thick = np.vstack((radii_I0_Thick, radii_I0_Thick_i))
    radii_I1_Thick = np.vstack((radii_I1_Thick, radii_I1_Thick_i))
    radii_I2_Thick = np.vstack((radii_I2_Thick, radii_I2_Thick_i))
    radii_FullAbsorption_Thick = np.vstack((radii_FullAbsorption_Thick, radii_FullAbsorption_Thick_i))

    # Total Flux Calcualtions
    # janksys_thin[i, 0] = ilp.total_jy(I0, brightparams["nu0"], brightparams["mass"]).value
    # janksys_thin[i, 1] = ilp.total_jy(I1, brightparams["nu0"], brightparams["mass"]).value
    # janksys_thin[i, 2] = ilp.total_jy(I2, brightparams["nu0"], brightparams["mass"]).value
    # janksys_thin[i, 3] = ilp.total_jy(I0 + I1 + I2, brightparams["nu0"], brightparams["mass"]).value

    janksys_thick[i, 0] = ilp.total_jy(I0_Absorb, brightparams["nu0"], brightparams["mass"]).value
    janksys_thick[i, 1] = ilp.total_jy(I1_Absorb, brightparams["nu0"], brightparams["mass"]).value
    janksys_thick[i, 2] = ilp.total_jy(I2_Absorb, brightparams["nu0"], brightparams["mass"]).value
    janksys_thick[i, 3] = ilp.total_jy(Absorbtion_Image, brightparams["nu0"], brightparams["mass"]).value

# Remove Row of Zeros
radii_Full_Thin = np.delete(radii_Full_Thin, 0, 0)
radii_I0_Thin = np.delete(radii_I0_Thin, 0, 0)
radii_I1_Thin = np.delete(radii_I1_Thin, 0, 0)
radii_I2_Thin = np.delete(radii_I2_Thin, 0, 0)

radii_FullAbsorption_Thick = np.delete(radii_FullAbsorption_Thick, 0, 0)
radii_I0_Thick = np.delete(radii_I0_Thick, 0, 0)
radii_I1_Thick = np.delete(radii_I1_Thick, 0, 0)
radii_I2_Thick = np.delete(radii_I2_Thick, 0, 0)

# Black Hole Inner Shadow Calc--------------------------
r_inner = np.load('r_inner_spin_{}_inc_{}.npy'.format(spin_case, i_case))
alphas_inner = np.load('alphas_inner_spin_{}_inc_{}.npy'.format(spin_case, i_case))
betas_inner = np.load('betas_inner_spin_{}_inc_{}.npy'.format(spin_case, i_case))

# Black Hole Outer Shadow Calc--------------------------
r_outer = np.load('r_outer_spin_{}_inc_{}.npy'.format(spin_case, i_case))
alphas_outer = np.load('alphas_outer_spin_{}_inc_{}.npy'.format(spin_case, i_case))
betas_outer = np.load('betas_outer_spin_{}_inc_{}.npy'.format(spin_case, i_case))

'''------------------------------------------------------------------------'''
jansky_variable_total = janksys_thick[:, 3]
jansky_variable_n0 = janksys_thick[:, 0]
jansky_variable_n1 = janksys_thick[:, 1]
jansky_variable_n2 = janksys_thick[:, 2]

ring_radii_n0_array = radii_I0_Thick
ring_radii_n1_array = radii_I1_Thick
ring_radii_n2_array = radii_I2_Thick

ring_radii_n0 = mean_radii_Thick[:, 0]
ring_radii_n1 = mean_radii_Thick[:, 1]
ring_radii_n2 = mean_radii_Thick[:, 2]

b = 0
# Figures--------------------------
for i in range(num_iterations):
    print('Creating Frame : ' + str(b))

    h5f = h5py.File(doth5_files[i], 'r')

    I0 = h5f['bghts0'][:]
    I1 = h5f['bghts1'][:]
    I2 = h5f['bghts2'][:]
    image = Absorbtion_Image = h5f['bghts_full_absorbtion'][:]

    h5f.close()

    # image-------------------------
    one_M = ilp.rg_func(brightparamsLEG[1] * u.g).to(u.m)
    M2uas = np.arctan(one_M.value / dBH) / muas_to_rad  # Mass to micro arcseconds

    # if b == 0:
    # 	vmax = np.max(image)

    vmax = np.nanmax(image)

    # fig = plt.subplots(1,1, figsize=[6,5],dpi=400)

    fig = plt.subplots(1, 2, figsize=[10, 5], dpi=400, width_ratios=[2, 1])

    xaxis = np.array(x_variable) / label[action[0], 1]

    ax = [None, None]

    ax[0] = plt.subplot(1, 2, 1)
    # im = ax[0].imshow(image,vmax=np.max(image)*1.2,origin="lower",cmap="afmhot",extent=[-lim0,lim0,-lim0,lim0])
    im = ax[0].imshow(image, origin="lower", cmap="afmhot", extent=[-lim0, lim0, -lim0, lim0],
                      norm=matplotlib.colors.PowerNorm(action[4], vmax=vmax))
    ax[0].set_xlim(-10, 10)  # units of M
    ax[0].set_ylim(-10, 10)

    ax[0].set_xlabel(r"$\alpha$" + " " + r"($\mu as$)")
    ax[0].set_ylabel(r"$\beta$" + " " + r"($\mu as$)")
    ax[0].text(-9, 8.5,
               label[action[0], 0] + str(round(x_variable[i] / label[action[0], 1], 2)) + ' ' + label[action[0], 2],
               fontsize=12, color="w")

    ax[0].set_xticks([-10, -7.5, -5, -2.5, 0, 2.5, 5, 7.5, 10], labels=[
        str('{:.3}'.format(-10 * M2uas)),
        str('{:.3}'.format(-7.5 * M2uas)),
        str('{:.3}'.format(-5 * M2uas)),
        str('{:.3}'.format(-2.5 * M2uas)),
        str('{:.3}'.format(0 * M2uas)),
        str('{:.3}'.format(2.5 * M2uas)),
        str('{:.3}'.format(5 * M2uas)),
        str('{:.3}'.format(7.5 * M2uas)),
        str('{:.3}'.format(10 * M2uas))
    ])

    ax[0].set_yticks([-10, -7.5, -5, -2.5, 0, 2.5, 5, 7.5, 10], labels=[
        str('{:.3}'.format(-10 * M2uas)),
        str('{:.3}'.format(-7.5 * M2uas)),
        str('{:.3}'.format(-5 * M2uas)),
        str('{:.3}'.format(-2.5 * M2uas)),
        str('{:.3}'.format(0 * M2uas)),
        str('{:.3}'.format(2.5 * M2uas)),
        str('{:.3}'.format(5 * M2uas)),
        str('{:.3}'.format(7.5 * M2uas)),
        str('{:.3}'.format(10 * M2uas))
    ])

    colorbar = plt.colorbar(im, fraction=0.046, pad=0.04, format='%.1e', ticks=[
        vmax * .8,
        vmax * .6,
        vmax * .4,
        vmax * .2,
        vmax * .05
    ],
                            label="Brightnes Temperature (K)"
                            )

    ax[1] = plt.subplot(1, 2, 2)
    # ax[1].axhline(r_inner, color='k', label='Blackhole Inner Shadow', linewidth=3)
    ax[1].axhline(r_outer, color='dimgrey', label='Blackhole Shadow', linewidth=3)
    ax[1].plot(xaxis, ring_radii_n0, '--', label='n=0', color='tab:red', linewidth=2)
    ax[1].plot(xaxis, ring_radii_n1, '--', label='n=1', color='tab:orange', linewidth=2)
    ax[1].plot(xaxis, ring_radii_n2, '--', label='n=2', color='tab:blue', linewidth=2)

    ax[1].set_xlabel(label[action[0], 0].replace('=', '') + ' (' + label[action[0], 2] + ')')
    ax[1].set_ylabel("Ring Radii ({})".format(R'$R_g$'))

    p5 = ax[1].scatter((xaxis)[i], ring_radii_n0[i], color='tab:red')
    p6 = ax[1].scatter((xaxis)[i], ring_radii_n1[i], color='tab:orange')
    p7 = ax[1].scatter((xaxis)[i], ring_radii_n2[i], color='tab:blue')

    ax[1].legend()
    ax[1].set_xlim(xaxis[0], xaxis[xaxis.size - 1])
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=.8, hspace=None)

    figname = 'Fig_{}.png'.format(b)

    # plt.colorbar(im)
    plt.savefig(figname, dpi=400, bbox_inches='tight')
    plt.close()
    names_to_delete2 += [figname]
    b = b + 1

movie_name = 'Present_BHMovie_var_{}_start_{}_stop_{}_steps_{}_a_{}_i_{}_nu_{}_mass_{}_scaleh_{}_thetab_{}_beta_{}_Rie_{}_Bchoi_{}_rb_{}_nth0_{}_te0_{}_pdens_{}_ptemp_{}.mp4'.format(
    cmd_args[action[0]],
    "{:.3e}".format(action[1]),  # start
    "{:.3e}".format(action[2]),  # stop
    "{:.3e}".format(action[3]),  # steps
    spin_case,  # Spin
    i_case,  # Observing angle
    "{:.1e}".format(brightparamsLEG[0]),  # Nu
    "{:.1e}".format(brightparamsLEG[1]),  # blkhole mass
    brightparamsLEG[2],  # scale height
    "{:.3e}".format(brightparamsLEG[3]),  # theta b
    "{:.1e}".format(brightparamsLEG[4]),  # beta
    "{:.1e}".format(brightparamsLEG[5]),  # Rie
    "{:.1e}".format(brightparamsLEG[6]),  # Bchoice
    "{:.1e}".format(brightparamsLEG[7]),  # rb
    "{:.1e}".format(brightparamsLEG[8]),  # nth0
    "{:.1e}".format(brightparamsLEG[9]),  # te0
    "{:.1e}".format(brightparamsLEG[10]),  # pdens
    "{:.1e}".format(brightparamsLEG[11])  # ptemp
)

if os.path.isfile('./' + movie_name):
    subprocess.run(['rm ' + './' + movie_name], shell=True)

speed = 8  # TODO: Check what units

subprocess.run(["ffmpeg -r " + str(
    speed) + " -i Fig_%d.png -vf 'pad=ceil(iw/2)*2:ceil(ih/2)*2' -vcodec libx264 -crf 10 -pix_fmt yuv420p " + movie_name],
               shell=True)

# with imageio.get_writer('BHMovie_var_{}_start_{}_stop_{}_steps_{}_a_{}_i_{}_nu_{}_mass_{}_scaleh_{}_thetab_{}_beta_{}_rb_{}_nth0_{}_te0_{}_pdens_{}_ptemp_{}.gif'.format(
#     cmd_args[action[0]], "{:.3e}".format(action[1]),"{:.3e}".format(action[2]),"{:.3e}".format(action[3]),spin_case,i_case,"{:.1e}".format(brightparams[0]),"{:.1e}".format(brightparams[1]), brightparams[2],
#     "{:.3e}".format(brightparams[3]), "{:.1e}".format(brightparams[4]),"{:.1e}".format(brightparams[5]), "{:.1e}".format(brightparams[6]),
#     "{:.1e}".format(brightparams[7]),"{:.1e}".format(brightparams[8]),"{:.1e}".format(brightparams[9])), mode='I') as writer:
#     for filename in names_to_delete2:
#         image = imageio.imread(filename)
#         writer.append_data(image)


for i in range(len(doth5_files)):
    subprocess.run(['rm ' + doth5_files[i]], shell=True)
    subprocess.run(['rm ' + names_to_delete2[i]], shell=True)
