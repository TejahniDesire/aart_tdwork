import sys
import subprocess

aartpath = '/home/td6241/repositories/aart'  # insert path to aart repo
sys.path.append(aartpath)
import image_tools
import kgeo

import classRunComputing

from aart_func import *
import params
from astropy import units as u

import runDataClass
import astroModels
import numpy as np
import movieMakerV2.movieMakerIntensity
from movieMakerV2 import movieMakerIntensity

'''Settings______________________________________'''

model_name = "ModelC22"
run = runDataClass.runData("run2",
                           astroModels.bp_run2,
                           ["p_temp", "p_mag"],
                           [(["a"], [str(.001)]), (["a"], [str(.5)]), (["a"], [str(15 / 16)])],
                           ["ModelA", "ModelB", "ModelC"],
                           )
run.setisNormalized(True)

action = run.getAction()

bigRun = classRunComputing.BigRuns(
    run.getRunName(),
    run.getBrightparams(),
    run.getBPVarNames(),
    run.getGeoGrid(),
    run.getGeoGridNames(),
    normalized_brightparams=run.getIsNormalized(),
)

sub_paths = bigRun.getSubPaths()
model = runDataClass.SingleModelData(sub_paths ,run.getRunName(),model_name)
brightparams = bigRun.getModelBrightParams(model_name)
num_iterations = int((action["stop"] - action["start"]) / action["step"])

data_path = sub_paths["intensityPath"] + model_name + "/clean/numpy/"

x_variable = np.load(data_path + "x_variable.npy")
janksys_thick = np.load(data_path + "janksys_thick.npy")
janksys_thin = np.load(data_path + "janksys_thin.npy")
mean_radii_Thin = np.load(data_path + "mean_radii_Thin.npy")
mean_radii_Thick = np.load(data_path + "mean_radii_Thick.npy")
radii_I0_Thin = np.load(data_path + "radii_I0_Thin.npy")
radii_I1_Thin = np.load(data_path + "radii_I1_Thin.npy")
radii_I2_Thin = np.load(data_path + "radii_I2_Thin.npy")
radii_Full_Thin = np.load(data_path + "radii_Full_Thin.npy")
radii_FullAbsorption_Thick = np.load(data_path + "radii_FullAbsorption_Thick.npy")
radii_I0_Thick = np.load(data_path + "radii_I0_Thick.npy")
radii_I1_Thick = np.load(data_path + "radii_I1_Thick.npy")
radii_I2_Thick = np.load(data_path + "radii_I2_Thick.npy")
theta = np.load(data_path + "theta.npy")
mean_optical_depth_I0 = np.load(data_path + "mean_optical_depth_I0.npy")
mean_optical_depth_I1 = np.load(data_path + "mean_optical_depth_I1.npy")
mean_optical_depth_I2 = np.load(data_path + "mean_optical_depth_I2.npy")

# Construct Shadows___________________________________________________________________
a = params.spin_case
inc = params.i_case * np.pi / 180  # inclination angle
rh = 1 + np.sqrt(1 - a ** 2)  # event horizon
# angles to sample
varphis = np.linspace(-180, 179, 360) * np.pi / 180

# generate inner shadow (n=0) curve with kgeo
data_inner = kgeo.equatorial_lensing.rho_of_req(a, inc, rh, mbar=0, varphis=varphis)
(_, rhos_inner, alphas_inner, betas_inner) = data_inner

r_inner = image_tools.curve_params(varphis, rhos_inner)

# generate outer shadow (n=inf) curve with kgeo
data_outer = kgeo.equatorial_lensing.rho_of_req(a, inc, rh, mbar=5, varphis=varphis)
(_, rhos_outer, alphas_outer, betas_outer) = data_outer

r_outer = image_tools.curve_params(varphis, rhos_outer)

one_M = ilp.rg_func(brightparams["mass"] * u.g).to(u.m)
M2uas = np.arctan(one_M.value / dBH) / muas_to_rad  # Mass to micro arcseconds

k = action["start"]
parent_model_path = sub_paths["intensityPath"] + model_name + "/"
current_model_file = parent_model_path + "clean/"

ring_radii_n0_array = radii_I0_Thick
ring_radii_n1_array = radii_I1_Thick
ring_radii_n2_array = radii_I2_Thick

ring_radii_n0 = mean_radii_Thick[:, 0]
ring_radii_n1 = mean_radii_Thick[:, 1]
ring_radii_n2 = mean_radii_Thick[:, 2]
names_to_delete2 = []  # images
for i in range(num_iterations):
    print('Creating Frame : ' + str(i))
    print(R"Observation frequency $\nu=$", k)
    # image-------------------------
    brightparams["nu0"] = k

    current_intensity_file = (current_model_file +
                              action["var"] + "_" + "{:.5e}".format(brightparams[action["var"]]))

    lim0 = 25

    print("Reading file: ", current_intensity_file)

    h5f = h5py.File(current_intensity_file, 'r')

    I0 = h5f['bghts0'][:]  # This implies I0 is 1 pass
    I1 = h5f['bghts1'][:]
    I2 = h5f['bghts2'][:]

    I2_Absorb = h5f['bghts2_absorbtion'][:]
    I1_Absorb = h5f['bghts1_absorbtion'][:]
    I0_Absorb = h5f['bghts0_absorbtion'][:]
    Absorbtion_Image = h5f['bghts_full_absorbtion'][:]

    h5f.close()

    image = Absorbtion_Image

    thick_radii = [radii_I0_Thick[i, :], radii_I1_Thick[i, :],
                   radii_I2_Thick[i, :], radii_FullAbsorption_Thick[i, :]]

    thick_alpha0 = thick_radii[0] * np.cos(theta)
    thick_beta0 = thick_radii[0] * np.sin(theta)
    thick_alpha1 = thick_radii[1] * np.cos(theta)
    thick_beta1 = thick_radii[1] * np.sin(theta)
    thick_alpha2 = thick_radii[2] * np.cos(theta)
    thick_beta2 = thick_radii[2] * np.sin(theta)
    thick_alpha_full = thick_radii[3] * np.cos(theta)
    thick_beta_full = thick_radii[3] * np.sin(theta)


    vmax = np.max(image)

    # fig = plt.subplots(1,1, figsize=[6,5],dpi=400)

    fig = plt.subplots(1, 2, figsize=[10, 5], dpi=400, width_ratios=[2, 1])


    ax = [None, None]

    ax[0] = plt.subplot(1, 2, 1)
    # im = ax[0].imshow(I0+I1+I2,vmax=np.max(I0+I1+I2)*1.2,origin="lower",cmap="afmhot",extent=[-lim0,lim0,-lim0,lim0])
    im = ax[0].imshow(image, origin="lower", cmap="afmhot", extent=[-lim0, lim0, -lim0, lim0],
                      norm=matplotlib.colors.PowerNorm(action[4], vmax=vmax))
    ax[0].set_xlim(-10, 10)  # units of M
    ax[0].set_ylim(-10, 10)

    ax[0].set_xlabel(r"$\alpha$" + " " + r"($\mu as$)")
    ax[0].set_ylabel(r"$\beta$" + " " + r"($\mu as$)")
    ax[0].text(-9, 8.5, astroModels.var_label[action["var"]]
             + str(round(x_variable[i] / astroModels.scale_label[action["var"]], 2))
             + ' ' + astroModels.units_label[action["var"]], fontsize=12, color="w")

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
    #
    # colorbar0 = fig.colorbar(im, fraction=0.046, pad=0.04, format=ticker.FuncFormatter(fmt), ticks=[
    #     vmax * .8,
    #     vmax * .6,
    #     vmax * .4,
    #     vmax * .2,
    #     vmax * .05
    # ],
    #                          ax=ax[0]
    #                          )

    ax[1] = plt.subplot(1, 2, 2)

    ax[1].plot(x_variable, mean_radii_Thick[:, 0], '-', label=R'n=0', color='tab:red', linewidth=3)
    ax[1].plot(x_variable, mean_radii_Thick[:, 1], ':', label=R'n=1', color='tab:orange', linewidth=3)
    ax[1].plot(x_variable, mean_radii_Thick[:, 2], '--', label=R'n=2', color='tab:blue', linewidth=3)
    ax[1].plot(x_variable, mean_radii_Thick[:, 3], '-.', label='Cumulative', color='tab:purple', linewidth=3)

    # ax[1].axhline(r_inner, color='k', label='Blackhole Inner Shadow', linewidth=3)
    ax[1].axhline(r_outer, color='dimgrey', label='Blackhole Shadow', linewidth=3)
    ax[1].plot(x_variable, ring_radii_n0, '--', label='n=0', color='tab:red', linewidth=2)
    ax[1].plot(x_variable, ring_radii_n1, '--', label='n=1', color='tab:orange', linewidth=2)
    ax[1].plot(x_variable, ring_radii_n2, '--', label='n=2', color='tab:blue', linewidth=2)

    # ax[1].set_xlabel(label[action[0], 0].replace('=', '') + ' (' + label[action[0], 2] + ')')
    ax[1].set_ylabel("Ring Radii ({})".format(R'$R_g$'))

    p5 = ax[1].scatter(x_variable[i], ring_radii_n0[i], color='tab:red')
    p6 = ax[1].scatter(x_variable[i], ring_radii_n1[i], color='tab:orange')
    p7 = ax[1].scatter(x_variable[i], ring_radii_n2[i], color='tab:blue')

    ax[1].legend()
    ax[1].set_xlim(x_variable[0], x_variable[x_variable.size - 1])
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=.8, hspace=None)

    figname = sub_paths["movie"] + 'Fig_{}.png'.format(i)

    # plt.colorbar(im)
    plt.savefig(figname, dpi=400, bbox_inches='tight')
    plt.close()
    names_to_delete2 += [figname]


movie_name =sub_paths["movie"] + str(model_name) + "_Movie.mp4"
if os.path.isfile('./' + movie_name):
    subprocess.run(['rm ' + './' + movie_name], shell=True)
speed = 8


subprocess.run(["ffmpeg -r " + str(
    speed) + " -i " + sub_paths["movie"] + "Fig_%d.png -vf 'pad=ceil(iw/2)*2:ceil(ih/2)*2' -vcodec libx264 -crf 10 -pix_fmt yuv420p " + movie_name],
               shell=True)






