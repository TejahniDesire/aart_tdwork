import sys
import EZPaths
from Movie_Maker import astroModels
import fileloading
import image_tools
import normalizingBrightparams

sys.path.append(EZPaths.aartPath)
import subprocess
from aart_func import *
import params
from params import *  # The file params.py contains all the relevant parameters for the simulations
import image_tools as tls
from astroModels import *
import numpy as np
from movieMakerV2 import intensityBlurr

speed = 8

line = "\n__________________________________________________\n"


def intensity_movie(action, sub_path, model: str, intent_grid_type, brightparams,frequency_list=None,funckey=funckeys):
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

    geo_model = model.split("_")[0]  # remove numbers from model
    lband = sub_path["GeoDoth5Path"] + geo_model + "Lensing" + ".h5"
    rtray = sub_path["GeoDoth5Path"] + geo_model + "RayTracing" + ".h5"
    MagAng = sub_path["GeoDoth5Path"] + geo_model + "MagAng" + ".h5"

    # File paths-------------------------
    parent_model_path = sub_path["intensityPath"] + model + "/"
    current_model_file = parent_model_path + "clean/"

    # total jy at 230GHz
    print("Finding Value at 230GHz...")
    thin_total_flux, thick_total_flux = normalizingBrightparams.totalIntensity230Point(lband, rtray,MagAng, brightparams,funckey,
                                                                                       already230=False)
    print("Value Found: ",thick_total_flux)
    intermodel_data = {
        "thin_total_flux": thin_total_flux,
        "thick_total_flux": thick_total_flux
    }

    # ----------------------------------------------------------------------------------------------------------------------
    # size = 1000  # array size for radii calcs
    # size = 500  # array size for radii calcs
    num_iterations = int((action["stop"] - action["start"]) / action["step"])

    """GRAPHS________________________________________________________________________________________________________"""
    # # Intensity images
    # I_thins = np.ndarray([num_iterations, 3])  # [I0, I1, I2]
    # I_thicks = np.ndarray([num_iterations, 3])  # [I0, I1, I2]
    # #

    done_list = None
    if frequency_list is not None:
        done_list = np.full(len(frequency_list),False)

    freq_xaxis = []
    for i in range(num_iterations):
        print(line)
        print(line)

        brightparams[action["var"]] = action["start"] + i * action["step"]
        current_freqeuncy = brightparams[action["var"]]
        freq_xaxis += [current_freqeuncy]
        do_image, done_list = fileloading.frequencyListAnalysis(frequency_list, done_list, current_freqeuncy)

        if do_image:
            print('Creating intensity.h5 for Model ' + model + ' number: ' + str(i))

            # create current intensity folder
            # Update varing parameter

            cmd = fileloading.createIntensityArgs(brightparams,lband,rtray,MagAng,funckey)
            subprocess.run([cmd], shell=True)

            # Read created file
            fnrays = fileloading.intensityNameNoUnits(brightparams,funckey)
            new_intensity_path = current_model_file + action["var"] + "_" + "{:.5e}".format(current_freqeuncy)
            subprocess.run(["mv " + fnrays + ' ' + new_intensity_path], shell=True)
        else:
            print('Skipping Intensity.h5 for Model ' + model + ' number: ' + str(i) + "...")
    intermodel_data["xaxis"] = freq_xaxis
    return intermodel_data


"""________________________________________________________________________________________________________"""


def profileAnalysis(action, sub_path, model: str, brightparams,frequency_list=None):
    parent_model_path = sub_path["intensityPath"] + model + "/"
    current_model_file = parent_model_path + "clean/"
    final_data_path = current_model_file + "numpy/"

    num_iterations = int((action["stop"] - action["start"]) / action["step"])

    done_list = None
    if frequency_list is not None:
        done_list = np.full(len(frequency_list), False)

    L = 0
    for i in range(num_iterations):
        print(line)
        print(line)
        print('Reading intensity.h5 for Model ' + model + ' number: ' + str(i))

        brightparams[action["var"]] = action["start"] + i * action["step"]
        current_freqeuncy = brightparams[action["var"]]

        do_image, done_list = fileloading.frequencyListAnalysis(frequency_list, done_list, current_freqeuncy)

        if do_image:
            print('Extracting profiles for Model ' + model + ' number: ' + str(L))

            intensity_path = current_model_file + action["var"] + "_" + "{:.5e}".format(current_freqeuncy)

            if frequency_list is not None:
                h5f = h5py.File(intensity_path, 'r')
                full_profiles0 = h5f['full_profiles0'][:]
                full_profiles1 = h5f['full_profiles1'][:]
                full_profiles2 = h5f['full_profiles2'][:]
                # full_profiles_unit = h5f['full_profiles_unit'][:]
                print("Frequency = " + str(current_freqeuncy / 1e9) +
                      " GHz for power law saving at desired frequency = " + str(frequency_list[L] / 1e9) + "GHz")

                np.save(final_data_path + "_full_profiles0_{}GHz".format("{:.5e}".format(frequency_list[L])),
                        full_profiles0)
                np.save(final_data_path + "_full_profiles1_{}GHz".format("{:.5e}".format(frequency_list[L])),
                        full_profiles1)
                np.save(final_data_path + "_full_profiles2_{}GHz".format("{:.5e}".format(frequency_list[L])),
                        full_profiles2)

                h5f.close()

            L += 1
        else:
            print('Skipping clean analysis for model ' + model + ' number: ' + str(i) + "...")


def imageAnalysis(action, sub_path, model: str, brightparams,average=True):
    parent_model_path = sub_path["intensityPath"] + model + "/"
    current_model_file = parent_model_path + "clean/"
    num_of_theta_points = image_tools.num_of_theta_points  # array size for radii calcs

    final_data_path = current_model_file + "numpy/"

    fileloading.creatSubDirectory(final_data_path, "final image path for {}".format(model), kill_policy=False)
    if not average:
        final_data_path += "FalseAvg_"

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
    mean_optical_depth_Total = np.zeros(num_iterations)

    for i in range(num_iterations):
        print(line)
        print(line)
        print('Reading intensity.h5 for Model ' + model + ' number: ' + str(i))

        brightparams[action["var"]] = action["start"] + i * action["step"]
        current_freqeuncy = brightparams[action["var"]]

        print('Analyzing intensity.h5 for Model ' + model + ' number: ' + str(i))
        x_variable[i] = current_freqeuncy

        intensity_path = current_model_file + action["var"] + "_" + "{:.5e}".format(current_freqeuncy)

        h5f = h5py.File(intensity_path, 'r')

        I0 = h5f['bghts0'][:]
        I1 = h5f['bghts1'][:]
        I2 = h5f['bghts2'][:]

        I0_Absorb = h5f['bghts0_absorbtion'][:]
        I1_Absorb = h5f['bghts1_absorbtion'][:]
        I2_Absorb = h5f['bghts2_absorbtion'][:]
        Absorbtion_Image = h5f['bghts_full_absorbtion'][:]

        tau2 = h5f['tau2'][()]
        tau1 = h5f['tau1'][()]
        tau0 = h5f['tau0'][()]
        tauTotal = h5f['tauTotal'][()]

        h5f.close()

        # Thin Radii Calcs----------------------------------------------------------------------------------------------

        radii_I0_Thin_i, theta = tls.radii_of_thetaV2(I0, params.dx0,average=average)
        radii_I1_Thin_i, theta = tls.radii_of_thetaV2(I1, params.dx0,average=average)
        radii_I2_Thin_i, theta = tls.radii_of_thetaV2(I2, params.dx0,average=average)
        radii_Full_Thick_i, theta = tls.radii_of_thetaV2(I0 + I1 + I2, params.dx0,average=average)

        # profs = scipy.ndimage.convolve1d(profs, np.ones(navg_ang), axis=0) / navg_ang

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
        radii_I0_Thick_i, theta = tls.radii_of_thetaV2(I0_Absorb, params.dx0,average=average)
        radii_I1_Thick_i, theta = tls.radii_of_thetaV2(I1_Absorb, params.dx0,average=average)
        radii_I2_Thick_i, theta = tls.radii_of_thetaV2(I2_Absorb, params.dx0,average=average)
        radii_FullAbsorption_Thick_i, theta = tls.radii_of_thetaV2(Absorbtion_Image, params.dx0,average=average)

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

        mean_optical_depth_I0[i] = tau0
        mean_optical_depth_I1[i] = tau1
        mean_optical_depth_I2[i] = tau2
        mean_optical_depth_Total[i] = tauTotal
        # (\Sum \tau * I) / (\Sum I),

    # Remove Row of Zeros
    radii_Full_Thin = np.delete(radii_Full_Thin, 0, 0)
    radii_I0_Thin = np.delete(radii_I0_Thin, 0, 0)
    radii_I1_Thin = np.delete(radii_I1_Thin, 0, 0)
    radii_I2_Thin = np.delete(radii_I2_Thin, 0, 0)

    radii_FullAbsorption_Thick = np.delete(radii_FullAbsorption_Thick, 0, 0)
    radii_I0_Thick = np.delete(radii_I0_Thick, 0, 0)
    radii_I1_Thick = np.delete(radii_I1_Thick, 0, 0)
    radii_I2_Thick = np.delete(radii_I2_Thick, 0, 0)

    # Saving Data-------------------------------------------------------------------------------------------------------
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
    np.save(final_data_path + "mean_optical_depth_I0", mean_optical_depth_I0)
    np.save(final_data_path + "mean_optical_depth_I1", mean_optical_depth_I1)
    np.save(final_data_path + "mean_optical_depth_I2", mean_optical_depth_I2)
    np.save(final_data_path + "mean_optical_depth_Total", mean_optical_depth_Total)


def blurr_intensity_movie(action, sub_path, model: str, intent_grid_type: int,
                          brightparams: dict, blurr_frequency_list, blur_kernal,funckey=funckeys):
    geo_model = model.split("_")[0]  # remove numbers from model
    lband = sub_path["GeoDoth5Path"] + geo_model + "Lensing" + ".h5"
    rtray = sub_path["GeoDoth5Path"] + geo_model + "RayTracing" + ".h5"
    MagAng = sub_path["GeoDoth5Path"] + geo_model + "MagAng" + ".h5"

    print("lband: '{}'".format(lband))
    print("rtray: '{}'".format(rtray))
    print("MagAng: '{}'".format(MagAng))

    # File paths-------------------------
    parent_model_path = sub_path["intensityPath"] + model + "/"
    current_model_file = parent_model_path + "blurr" + str(blur_kernal) + "/"
    clean_model_file = parent_model_path + "clean/"

    fileloading.creatSubDirectory(current_model_file,
                                  "for {} blurr intensities".format(model), kill_policy=True)

    # total jy at 230GHz
    print("Calculating total flux at 230GHz")
    thin_total_230flux, thick_total_230flux = normalizingBrightparams.totalIntensity230Point(
        lband, rtray, MagAng, brightparams,funckey, already230=False, blurr_policy=True,blur_kernal=blur_kernal)

    intermodel_data = {
        "thin_total_flux": thin_total_230flux,
        "thick_total_flux": thick_total_230flux
    }

    # ----------------------------------------------------------------------------------------------------------------------
    # size = 1000  # array size for radii calcs
    # size = 500  # array size for radii calcs
    num_iterations = int((action["stop"] - action["start"]) / action["step"])

    done_list = None
    if blurr_frequency_list is not None:
        done_list = np.full(len(blurr_frequency_list),False)

    for i in range(num_iterations):
        print(line)
        print(line)
        brightparams[action["var"]] = action["start"] + i * action["step"]
        current_freqeuncy = brightparams[action["var"]]

        do_blurr, done_list = fileloading.frequencyListAnalysis(blurr_frequency_list, done_list, current_freqeuncy)

        if do_blurr:
            print('Reading intensity.h5 for Model ' + model + ' number: ' + str(i))
            clean_intensity_path = clean_model_file + action["var"] + "_" + "{:.5e}".format(current_freqeuncy)

            # Read File
            h5f = h5py.File(clean_intensity_path, 'r')

            I0 = h5f['bghts0'][:]
            I1 = h5f['bghts1'][:]
            I2 = h5f['bghts2'][:]

            absorb_image = h5f['bghts_full_absorbtion'][:]
            thin_image = I0 + I1 + I2

            h5f.close()

            # Blurring ______________________________________
            thin_blurr_image, thick_blurr_image = intensityBlurr.blurrIntensity(
                brightparams["mass"], thin_image, absorb_image,blur_kernal)
            # ______________________________________

            blurr_intensity_path = (current_model_file +
                                    action["var"] + "_blurr_" + "{:.5e}".format(current_freqeuncy))
            h5f = h5py.File(blurr_intensity_path, 'w')
            h5f.create_dataset('thin_blurr_image', data=thin_blurr_image)
            h5f.create_dataset("thick_blurr_image", data=thick_blurr_image)
            h5f.close()
            print("File ", blurr_intensity_path, " created.")
        else:
            print("Freuquency {} marked for skipping...".format("{:.5e}".format(current_freqeuncy)))

    return intermodel_data


def blurrImageAnalysis(action, sub_path, model: str, brightparams,blurr_frequency_list,blur_kernal):
    parent_model_path = sub_path["intensityPath"] + model + "/"
    current_model_file = parent_model_path + "blurr" + str(blur_kernal) + "/"
    num_of_theta_points = image_tools.num_of_theta_points  # array size for radii calcs
    num_iterations = int((action["stop"] - action["start"]) / action["step"])

    if blurr_frequency_list is None:
        len_of_freq_list = num_iterations
    else:
        len_of_freq_list = len(blurr_frequency_list)

    """GRAPH Variables________________________________________________________________________________________________________"""
    x_variable = np.zeros(len_of_freq_list)  # counter for independant variable

    # flux_________________________________
    janksys_thick = np.ndarray([len_of_freq_list, 1])  # [FullImage]
    janksys_thin = np.ndarray([len_of_freq_list, 1])  # [total]

    # Average Radius over all Theta_________________________________
    mean_radii_Thin = np.ndarray([len_of_freq_list, 1])  # [FullImage]
    mean_radii_Thick = np.ndarray([len_of_freq_list, 1])  # [FullImage]

    # Radius as a function of theta_________________________________
    radii_cumulative_Thin = np.zeros(num_of_theta_points)

    radii_cumulative_Thick = np.zeros(num_of_theta_points)

    """________________________________________________________________________________________________________"""

    done_list = None
    if blurr_frequency_list is not None:
        done_list = np.full(len(blurr_frequency_list),False)

    L = 0
    for i in range(num_iterations):
        print('Reading intensity.h5 for Model ' + model + ' number: ' + str(i))
        print(line)
        print(line)

        brightparams[action["var"]] = action["start"] + i * action["step"]
        current_freqeuncy = brightparams[action["var"]]

        do_blurr, done_list = fileloading.frequencyListAnalysis(blurr_frequency_list, done_list, current_freqeuncy)

        if do_blurr:
            print('Reading blurred intensity.h5 for Model ' + model + ' number: ' + str(L))
            x_variable[L] = current_freqeuncy
            intensity_path = current_model_file + action["var"] + "_blurr_" + "{:.5e}".format(current_freqeuncy)

            h5f = h5py.File(intensity_path, 'r')

            Thin_Image = h5f['thin_blurr_image'][:]
            Absorbtion_Image = h5f["thick_blurr_image"][:]

            h5f.close()

            # Thin Radii Calcs------------------------------------------------------------------------------------------

            radii_Thin_i, theta = tls.radii_of_thetaV2(Thin_Image, params.dx0)

            r0_thin = tls.curve_params(theta, radii_Thin_i)

            mean_radii_Thin[L, 0] = r0_thin

            radii_cumulative_Thin = np.vstack((radii_cumulative_Thin, radii_Thin_i))

            # Thick Radii Calcs-----------------------------------------------------------------------------------------
            radii_FullAbsorption_Thick_i, theta = tls.radii_of_thetaV2(Absorbtion_Image, params.dx0)

            full_thick = tls.curve_params(theta, radii_FullAbsorption_Thick_i)

            mean_radii_Thick[L, 0] = full_thick

            radii_cumulative_Thick = np.vstack((radii_cumulative_Thick, radii_FullAbsorption_Thick_i))

            # Total Flux Calcualtions
            janksys_thin[L, 0] = ilp.total_jy(Thin_Image, brightparams["nu0"], brightparams["mass"]).value

            janksys_thick[L, 0] = ilp.total_jy(Absorbtion_Image, brightparams["nu0"], brightparams["mass"]).value

            L += 1

        else:
            print("Freuquency {} marked for skipping...".format("{:.5e}".format(current_freqeuncy)))

    final_data_path = current_model_file + "numpy/"

    fileloading.creatSubDirectory(final_data_path, "final image path for {}".format(model), kill_policy=False)

    # Remove Row of Zeros
    radii_cumulative_Thin = np.delete(radii_cumulative_Thin, 0, 0)
    radii_cumulative_Thick = np.delete(radii_cumulative_Thick, 0, 0)

    # Saving Data--------------------------------------------------------------------------------------------------------
    np.save(final_data_path + "blurr_x_variable", x_variable)
    np.save(final_data_path + "blurr_janksys_thick", janksys_thick)
    np.save(final_data_path + "blurr_janksys_thin", janksys_thin)
    np.save(final_data_path + "blurr_mean_radii_Thin", mean_radii_Thin)
    np.save(final_data_path + "blurr_mean_radii_Thick", mean_radii_Thick)
    np.save(final_data_path + "blurr_radii_I0_Thin", radii_cumulative_Thin)
    np.save(final_data_path + "blurr_radii_FullAbsorption_Thick", radii_cumulative_Thick)
    np.save(final_data_path + "blurr_theta", theta)
