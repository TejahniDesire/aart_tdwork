import subprocess

import kgeo
import numpy as np
from matplotlib import ticker


import EZPaths
import os

import astroPloting
import bigRunComputing
import image_tools
from aart_func import *
from image_tools import curve_params
from params import *
import importlib
import params
import astroModels
import fileloading
from movieMakerV2 import movieMakerIntensity
import normalizingBrightparams
from astropy import units as u
import re


def numOfModel(model:str,all_full_model_names):
    j = 0
    for i in range(len(all_full_model_names)):
        if all_full_model_names[i] == model:
            print(all_full_model_names[i] + " Identified as " + model + " for number " + str(j))
            return j
        j += 1


def playModel(sub_path,save_path, run,action, model:str, intent_grid_type=2):

    print("Running " + model)
    print(bigRunComputing.line)
    print(bigRunComputing.line)
    print("Initializing Graph Creation")
    all_full_model_names = np.load(sub_path["meta"] + "AllModelsList.npy")
    all_brightparams = np.load(sub_path["meta"] + "AllBrightParamsList.npy", allow_pickle=True)
    thin_total_flux = np.load(sub_path["meta"] + "thin_total_flux.npy")
    thick_total_flux = np.load(sub_path["meta"] + "thick_total_flux.npy")

    j = numOfModel(model,all_full_model_names)

    current_geo_model = model[0:len(model) - intent_grid_type]
    fileloading.loadGeoModel(current_geo_model, run)
    lband = sub_path["GeoDoth5Path"] + current_geo_model + "Lensing" + ".h5"
    rtray = sub_path["GeoDoth5Path"] + current_geo_model + "RayTracing" + ".h5"

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
    # ___________________________________________________________________


    '''Data Readind----------------------------------'''
    data_path = sub_path["intensityPath"] + model + "/" + "numpy/"

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

    num_of_intensity_points = janksys_thin[:, 0].shape[0]
    print("Number of Intensity Points: ", num_of_intensity_points)

    k = action["start"]
    print("Constructing Full images for " + model)
    for i in range(num_of_intensity_points):
        brightparams = all_brightparams[j]
        brightparams["nu0"] = k
        print("Full image production for intensity frame: ", i)
        print(R"Observation frequency $\nu=$", k)

        current_intensity_file = (sub_path["intensityPath"] + model + "/" + action["var"]
                                  + "_" + "{:.5e}".format(brightparams[action["var"]]))

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

        thin_intensity = [I0, I1, I2, I0 + I1 + I2]
        thick_intensity = [I0_Absorb, I1_Absorb, I2_Absorb, Absorbtion_Image]
        rmax = I0.shape[0] * .4

        dim = [15, 7]
        '''IntensityVSRadiiType1________________________________________________________________'''
        fig, dum = plt.subplots(2, 2, figsize=dim, dpi=400)
        ax0 = plt.subplot(2, 2, 1)
        ax1 = plt.subplot(2, 2, 2)
        ax2 = plt.subplot(2, 2, 3)
        ax3 = plt.subplot(2, 2, 4)

        astroPloting.IntensityVSRadiiType1(fig, ax0, ax1,ax2,ax3,params.limits,thin_intensity,thick_intensity,rmax)

        pltname = (save_path['intVRad'] + 'IntVRad_' + str(i) + "_Nu_"
                   + str(round(x_variable[i] / astroModels.scale_label[action["var"]], 2)) + ".jpeg")
        plt.savefig(pltname, bbox_inches='tight')
        print("Image '{}' Created".format(pltname))
        plt.close()

        '''IntensityVSRadiiType2________________________________________________________________'''
        fig, dum = plt.subplots(1, 2, figsize=dim, dpi=400)
        ax0 = plt.subplot(1, 2, 1)
        ax1 = plt.subplot(1, 2, 2)

        astroPloting.IntensityVSRadiiType2(fig, ax0, ax1,params.limits,thin_intensity,rmax)

        pltname = (save_path['intVRad2'] + 'IntVRad_' + str(i) + "_Nu_"
                   + str(round(x_variable[i] / astroModels.scale_label[action["var"]], 2)) + ".jpeg")
        plt.savefig(pltname, bbox_inches='tight')
        print("Image '{}' Created".format(pltname))
        plt.close()

        '''RadVSVarphiType2________________________________________________________________'''
        fig, dum = plt.subplots(1, 2, figsize=dim, dpi=400)
        ax0 = plt.subplot(1, 2, 1)
        ax1 = plt.subplot(1, 2, 2)

        astroPloting.radiiVSVarphi(fig, ax0, ax1,params.limits,thin_intensity)

        pltname = (save_path['radVVarphi'] + 'radVVarphu_' + str(i) + "_Nu_"
                   + str(round(x_variable[i] / astroModels.scale_label[action["var"]], 2)) + ".jpeg")
        plt.savefig(pltname, bbox_inches='tight')
        print("Image '{}' Created".format(pltname))
        plt.close()

        k += action['step']



