import os.path
import sys
import EZPaths
import astroModels
import bigRunComputing
import fileloading
import image_tools

sys.path.append(EZPaths.aartPath)
import subprocess
import scipy
from matplotlib import ticker
from aart_func import *
import params
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
    else:
        subprocess.run(["rm -r " + current_model_file], shell=True)
        print("Subfolder {} already exist, deleting".format(current_model_file))
        subprocess.run(["mkdir " + current_model_file], shell=True)

    print("Subfolder for {} intensities created ({})".format(model,current_model_file))

    # ----------------------------------------------------------------------------------------------------------------------
    # size = 1000  # array size for radii calcs
    # size = 500  # array size for radii calcs
    size = image_tools.size  # array size for radii calcs
    num_iterations = int((action["stop"] - action["start"]) / action["step"])

    """GRAPHS________________________________________________________________________________________________________"""
    x_variable = np.zeros(num_iterations)  # counter for independant variable

    # flux_________________________________
    janksys_thick = np.ndarray([num_iterations, 4])  # [I0, I1, I2, FullImage]
    janksys_thin = np.ndarray([num_iterations, 4])  # [I0, I1, I2, total]

    # Average Radius over all Theta_________________________________
    mean_radii_Thin = np.ndarray([num_iterations, 4])  # [I0, I1, I2, FullImage]
    mean_radii_Thick = np.ndarray([num_iterations, 4])  # [I0, I1, I2, FullImage]

    # First layer of Radius as a function of theta_________________________________
    radii_Full_Thin = np.zeros(size)
    radii_I0_Thin = np.zeros(size)
    radii_I1_Thin = np.zeros(size)
    radii_I2_Thin = np.zeros(size)

    radii_FullAbsorption_Thick = np.zeros(size)
    radii_I0_Thick = np.zeros(size)
    radii_I1_Thick = np.zeros(size)
    radii_I2_Thick = np.zeros(size)

    # Optical Depth_________________________________
    mean_optical_depth_I0 = np.zeros(num_iterations)
    mean_optical_depth_I1 = np.zeros(num_iterations)
    mean_optical_depth_I2 = np.zeros(num_iterations)

    # total jy at 230GHz

    thin_total_flux, thick_total_flux = totalIntensity230Point(sub_path, model, intent_grid_type, brightparams)

    intermodel_data = {
        "thin_total_flux": thin_total_flux,
        "thick_total_flux":thick_total_flux
    }

    # # Intensity images
    # I_thins = np.ndarray([num_iterations, 3])  # [I0, I1, I2]
    # I_thicks = np.ndarray([num_iterations, 3])  # [I0, I1, I2]
    # #
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

        fnrays = fileloading.intensityNameNoUnits(brightparams, astroModels.funckeys)
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

        tau2 = h5f['tau2'][:]
        tau1 = h5f['tau1'][:]
        tau0 = h5f['tau0'][:]
        # full_profiles0 = h5f['full_profiles0'][:]
        # full_profiles1 = h5f['full_profiles1'][:]
        # full_profiles2 = h5f['full_profiles2'][:]
        # full_profiles_unit = h5f['full_profiles_unit'][:]


        # Thin Radii Calcs----------------------------------------------------------------------------------------------

        radii_I0_Thin_i, theta = tls.radii_of_thetaV2(I0, params.dx0)
        radii_I1_Thin_i, theta = tls.radii_of_thetaV2(I1, params.dx0)
        radii_I2_Thin_i, theta = tls.radii_of_thetaV2(I2, params.dx0)
        radii_Full_Thick_i, theta = tls.radii_of_thetaV2(I0 + I1 + I2, params.dx0)

        r0_thin = tls.curve_params(theta, radii_I0_Thin_i)
        r1_thin = tls.curve_params(theta, radii_I1_Thin_i)
        r2_thin = tls.curve_params(theta, radii_I2_Thin_i)
        full_thin = tls.curve_params(theta, radii_Full_Thick_i)

        mean_radii_Thin[i, 0] = r0_thin
        mean_radii_Thin[i, 1] = r1_thin
        mean_radii_Thin[i, 2] = r2_thin
        mean_radii_Thin[i, 3] = full_thin

        radii_I0_Thin = np.vstack((radii_I0_Thin, radii_I0_Thin_i))
        radii_I1_Thin = np.vstack((radii_I1_Thin, radii_I1_Thin_i))
        radii_I2_Thin = np.vstack((radii_I2_Thin, radii_I2_Thin_i))
        radii_Full_Thin = np.vstack((radii_Full_Thin, radii_Full_Thick_i))

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
        janksys_thin[i, 0] = ilp.total_jy(I0, brightparams["nu0"], brightparams["mass"]).value
        janksys_thin[i, 1] = ilp.total_jy(I1, brightparams["nu0"], brightparams["mass"]).value
        janksys_thin[i, 2] = ilp.total_jy(I2, brightparams["nu0"], brightparams["mass"]).value
        janksys_thin[i, 3] = ilp.total_jy(I0 + I1 + I2, brightparams["nu0"], brightparams["mass"]).value

        janksys_thick[i, 0] = ilp.total_jy(I0_Absorb, brightparams["nu0"], brightparams["mass"]).value
        janksys_thick[i, 1] = ilp.total_jy(I1_Absorb, brightparams["nu0"], brightparams["mass"]).value
        janksys_thick[i, 2] = ilp.total_jy(I2_Absorb, brightparams["nu0"], brightparams["mass"]).value
        janksys_thick[i, 3] = ilp.total_jy(Absorbtion_Image, brightparams["nu0"], brightparams["mass"]).value


        # Optical Depth-------------------------------------------------------------------------------------------------

        mean_optical_depth_I0[i] = np.sum(tau0 * I0) / np.sum(I0)
        mean_optical_depth_I1[i] = np.sum(tau0 * I1) / np.sum(I1)
        mean_optical_depth_I2[i] = np.sum(tau0 * I2) / np.sum(I2)
        # (\Sum \tau * I) / (\Sum I),

        h5f.close()



    final_data_path = current_model_file + "numpy/"

    if not os.path.isdir(final_data_path):
        subprocess.run(["mkdir " + final_data_path], shell=True)

    # Remove Row of Zeros
    radii_Full_Thin = np.delete(radii_Full_Thin, 0, 0)
    radii_I0_Thin = np.delete(radii_I0_Thin, 0, 0)
    radii_I1_Thin = np.delete(radii_I1_Thin, 0, 0)
    radii_I2_Thin = np.delete(radii_I2_Thin, 0, 0)

    radii_FullAbsorption_Thick = np.delete(radii_FullAbsorption_Thick, 0, 0)
    radii_I0_Thick = np.delete(radii_I0_Thick, 0, 0)
    radii_I1_Thick = np.delete(radii_I1_Thick, 0, 0)
    radii_I2_Thick = np.delete(radii_I2_Thick, 0, 0)

    # Saving Data--------------------------------------------------------------------------------------------------------
    np.save(final_data_path + "x_variable", x_variable)
    np.save(final_data_path + "janksys_thick", janksys_thick)
    np.save(final_data_path + "janksys_thin", janksys_thin)
    np.save(final_data_path + "mean_radii_Thin", mean_radii_Thin)
    np.save(final_data_path + "mean_radii_Thick", mean_radii_Thick)
    np.save(final_data_path + "radii_I0_Thin", radii_I0_Thin)
    np.save(final_data_path + "radii_I1_Thin", radii_I1_Thin)
    np.save(final_data_path + "radii_I2_Thin", radii_I2_Thin)
    np.save(final_data_path + "radii_Full_Thin", radii_Full_Thin)
    np.save(final_data_path + "radii_FullAbsorption_Thick", radii_FullAbsorption_Thick)
    np.save(final_data_path + "radii_I0_Thick", radii_I0_Thick)
    np.save(final_data_path + "radii_I1_Thick", radii_I1_Thick)
    np.save(final_data_path + "radii_I2_Thick", radii_I2_Thick)
    np.save(final_data_path + "theta", theta)
    np.save(final_data_path + "mean_optical_depth_I0",mean_optical_depth_I0)
    np.save(final_data_path + "mean_optical_depth_I1", mean_optical_depth_I1)
    np.save(final_data_path + "mean_optical_depth_I2", mean_optical_depth_I2)

    return intermodel_data


def totalIntensity230Point(sub_path, model:str, intent_grid_type, brightparams:dict):
    bp = brightparams.copy()
    bp["nu0"] = 230e9

    geo_model = model[0:len(model)-intent_grid_type]  # remove numbers from model
    lband = sub_path["GeoDoth5Path"] + geo_model + "Lensing" + ".h5"
    rtray = sub_path["GeoDoth5Path"] + geo_model + "RayTracing" + ".h5"

    args = bigRunComputing.createIntensityArgs(bp)
    args += "--lband " + lband + " --rtray " + rtray

    subprocess.run(['python3 ' + EZPaths.aartPath + '/radialintensity.py' + args], shell=True)

    # Read created file

    fnrays = fileloading.intensityNameNoUnits(brightparams, astroModels.funckeys)
    h5f = h5py.File(fnrays, 'r')

    I0 = h5f['bghts0'][:]
    I1 = h5f['bghts1'][:]
    I2 = h5f['bghts2'][:]

    # I0_Absorb = h5f['bghts0_absorbtion'][:]
    # I1_Absorb = h5f['bghts1_absorbtion'][:]
    # I2_Absorb = h5f['bghts2_absorbtion'][:]
    Absorbtion_Image = h5f['bghts_full_absorbtion'][:]
    #
    # tau2 = h5f['tau2'][:]
    # tau1 = h5f['tau1'][:]
    # tau0 = h5f['tau0'][:]

    h5f.close()

    subprocess.run(["rm " + fnrays], shell=True)

    thin_total_flux = ilp.total_jy(I0 + I1 + I2,230e9,bp["mass"])
    thick_total_flux = ilp.total_jy(Absorbtion_Image,230e9,bp["mass"])

    return thin_total_flux,thick_total_flux



