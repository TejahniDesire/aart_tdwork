import os.path
import sys
import EZPaths
import bigRunComputing
import fileloading
import image_tools

sys.path.append(EZPaths.aartPath)
import subprocess
import scipy
from matplotlib import ticker
from aart_func import *
from params import *  # The file params.py contains all the relevant parameters for the simulations
from astropy import units as u
import image_tools as tls
import numpy as np


speed = 8

line = "\n__________________________________________________\n"


def intensity_movie(action,sub_path, model:str, intent_grid_type, brightparams):
    """
        sub_paths = {
        "GeoDoth5Path"
        "intensityPath"
        "fluxPath"
        "radPath"
        "imagePath"

         action = {
         "var":
         "start":
         "stop":
         "step":
    """

    geo_model = model[0:len(model)-intent_grid_type] # remove numbers from model
    lband = sub_path["GeoDoth5Path"] + geo_model + "Lensing" + ".h5"
    rtray = sub_path["GeoDoth5Path"] + geo_model + "RayTracing" + ".h5"
    # '''Reading of the lensing bands----------------------------------'''
    # fnbands = lband
    #
    # '''Reading Analytical Ray-tracing----------------------------------'''
    # fnrays = rtray


    # '''Computing images----------------------------------'''

    brightparams_label = dict(brightparams)
    # brightparam_items = brightparams.items()
    # brightparams_init = np.array(list(brightparam_items))[:,1]  # Convert object to a list

    # File paths-------------------------
    current_model_file = sub_path["intensityPath"] + model + "/"

    if not os.path.isdir(current_model_file):
        subprocess.run(["mkdir " + current_model_file], shell=True)

    print("Subfolder for {} intensities created ({})".format(model,current_model_file))

    # ----------------------------------------------------------------------------------------------------------------------
    # size = 1000  # array size for radii calcs
    # size = 500  # array size for radii calcs
    size = image_tools.size  # array size for radii calcs
    num_iterations = int((action["stop"] - action["start"]) / action["step"])

    # Graphs
    x_variable = np.zeros(num_iterations)  # counter for independant variable

    janksys_thick = np.ndarray([num_iterations, 4])  # [I0, I1, I2, FullImage]
    janksys_thin = np.ndarray([num_iterations, 4])  # [I0, I1, I2, total]

    mean_radii_Thin = np.ndarray([num_iterations, 3])  # [I0, I1, I2]
    mean_radii_Thick = np.ndarray([num_iterations, 4])  # [I0, I1, I2, FullImage]

    radii_I0_Thin = np.zeros(size)
    radii_I1_Thin = np.zeros(size)
    radii_I2_Thin = np.zeros(size)

    radii_FullAbsorption_Thick = np.zeros(size)
    radii_I0_Thick = np.zeros(size)
    radii_I1_Thick = np.zeros(size)
    radii_I2_Thick = np.zeros(size)

    # # Intensity images
    # I_thins = np.ndarray([num_iterations, 3])  # [I0, I1, I2]
    # I_thicks = np.ndarray([num_iterations, 3])  # [I0, I1, I2]
    # #
    funckeys = {
        "emodelkey": 0,  # emodelkey Emission Model choice, 0 = thermal ultrarelativistic, 1 = power law
        "bkey": 2,  # bkey
        "nnoisykey": 0,  # nnoisykey Inoisy density. 0 = no noise, 1 = noise
        "tnoisykey": 0,  # tnoisykey Inoisy temperature
        "bnoisykey": 0  # bnoisykey Inoisy magnetic field
    }

    for i in range(num_iterations):
        print(line)
        print(line)
        print('Creating intensity.h5 for Model ' + model + ' number: ' + str(i))


        current_int_file = sub_path["intensityPath"]

        # create current intensity folder
        # Update varing parameter
        brightparams[action["var"]] = action["start"] + i * action["step"]
        x_variable[i] = brightparams[action["var"]]

        args = bigRunComputing.createIntensityArgs(brightparams)
        args += "--lband " + lband + " --rtray " + rtray


        subprocess.run(['python3 ' + EZPaths.aartPath + '/radialintensity.py' + args], shell=True)

        # Read created file

        fnrays = fileloading.intensityNameNoUnits(brightparams, funckeys)
        new_intensity_path = current_model_file + action["var"] + "_" + "{:.5e}".format(brightparams[action["var"]])
        subprocess.run(["mv " + fnrays + ' ' + new_intensity_path], shell=True)


        h5f = h5py.File(new_intensity_path, 'r')

        I0 = h5f['bghts0'][:]
        I1 = h5f['bghts1'][:]
        I2 = h5f['bghts2'][:]

        I0_Absorb = h5f['bghts0_absorbtion'][:]
        I1_Absorb = h5f['bghts1_absorbtion'][:]
        I2_Absorb = h5f['bghts2_absorbtion'][:]
        Absorbtion_Image = h5f['bghts_full_absorbtion'][:]

        # tau2 = h5f['tau2'][:]
        # tau1 = h5f['tau1'][:]
        # tau0 = h5f['tau0'][:]
        # full_profiles0 = h5f['full_profiles0'][:]
        # full_profiles1 = h5f['full_profiles1'][:]
        # full_profiles2 = h5f['full_profiles2'][:]
        # full_profiles_unit = h5f['full_profiles_unit'][:]


        # Thin Radii Calcs--------------------------------------------------------------------------------------------------

        radii_I0_Thin_i, theta = tls.radii_of_theta(I0, size)
        radii_I1_Thin_i, theta = tls.radii_of_theta(I1, size)
        radii_I2_Thin_i, theta = tls.radii_of_theta(I2, size)

        r0_thin = tls.curve_params(theta, radii_I0_Thin_i)
        r1_thin = tls.curve_params(theta, radii_I1_Thin_i)
        r2_thin = tls.curve_params(theta, radii_I2_Thin_i)

        mean_radii_Thin[i, 0] = r0_thin
        mean_radii_Thin[i, 1] = r1_thin
        mean_radii_Thin[i, 2] = r2_thin

        radii_I0_Thin = np.vstack((radii_I0_Thin, radii_I0_Thin_i))
        radii_I1_Thin = np.vstack((radii_I1_Thin, radii_I1_Thin_i))
        radii_I2_Thin = np.vstack((radii_I2_Thin, radii_I2_Thin_i))

        # Thick Radii Calcs-------------------------------------------------------------------------------------------------
        radii_I0_Thick_i, theta = tls.radii_of_theta(I0_Absorb, size)
        radii_I1_Thick_i, theta = tls.radii_of_theta(I1_Absorb, size)
        radii_I2_Thick_i, theta = tls.radii_of_theta(I2_Absorb, size)
        radii_FullAbsorption_Thick_i, theta = tls.radii_of_theta(Absorbtion_Image, size)

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
        janksys_thin[i, 0] = ilp.total_jy(I0, brightparams["nu0"], brightparams["mass"]).value
        janksys_thin[i, 1] = ilp.total_jy(I1, brightparams["nu0"], brightparams["mass"]).value
        janksys_thin[i, 2] = ilp.total_jy(I2, brightparams["nu0"], brightparams["mass"]).value
        janksys_thin[i, 3] = ilp.total_jy(I0 + I1 + I2, brightparams["nu0"], brightparams["mass"]).value

        janksys_thick[i, 0] = ilp.total_jy(I0_Absorb, brightparams["nu0"], brightparams["mass"]).value
        janksys_thick[i, 1] = ilp.total_jy(I1_Absorb, brightparams["nu0"], brightparams["mass"]).value
        janksys_thick[i, 2] = ilp.total_jy(I2_Absorb, brightparams["nu0"], brightparams["mass"]).value
        janksys_thick[i, 3] = ilp.total_jy(Absorbtion_Image, brightparams["nu0"], brightparams["mass"]).value

        h5f.close()

    final_data_path = current_model_file + "numpy/"

    if not os.path.isdir(final_data_path):
        subprocess.run(["mkdir " + final_data_path], shell=True)

    # Saving Data--------------------------------------------------------------------------------------------------------
    np.save(final_data_path + "x_variable", x_variable)
    np.save(final_data_path + "janksys_thick", janksys_thick)
    np.save(final_data_path + "janksys_thin", janksys_thin)
    np.save(final_data_path + "mean_radii_Thin", mean_radii_Thin)
    np.save(final_data_path + "mean_radii_Thick", mean_radii_Thick)
    np.save(final_data_path + "radii_I0_Thin", radii_I0_Thin)
    np.save(final_data_path + "radii_I1_Thin", radii_I1_Thin)
    np.save(final_data_path + "radii_I2_Thin", radii_I2_Thin)
    np.save(final_data_path + "radii_FullAbsorption_Thick", radii_FullAbsorption_Thick)
    np.save(final_data_path + "radii_I0_Thick", radii_I0_Thick)
    np.save(final_data_path + "radii_I1_Thick", radii_I1_Thick)
    np.save(final_data_path + "radii_I2_Thick", radii_I2_Thick)
    np.save(final_data_path + "theta", theta)