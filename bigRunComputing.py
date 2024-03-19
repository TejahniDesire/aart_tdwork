import subprocess

import kgeo
import numpy as np
from matplotlib import ticker


import EZPaths
import os

import astroPloting
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

kw_action = {
    "var": "nu0",
    "start": 3.00e+10,
    "stop": 6.00e+10,
    "step": 1.00e+10
}
funckeys = astroModels.funckeys


def createIntensityArgs(brightparams):
    args = ' '
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
        "b_0": '--b0 ',
        "p_dens": '--pdens ',
        "p_temp": '--ptemp ',
        "p_mag": '--pmag ',
        "nscale": '--nscale ',
    }
    cmd2_args = {
        "emodelkey": '--emodelkey ',
        "bkey": '--bkey ',
        "nnoisykey": '--nnoisykey ',
        "tnoisykey": '--tnoisykey ',
        "bnoisykey": '--bnoisykey ',
    }


    # brightparams = fpp.bp_steeperT
    # funckeys = fpp.fk_fiducial

    for arg in cmd1_args:
        args = args + cmd1_args[arg] + str(brightparams[arg]) + ' '

    for arg in cmd2_args:
        args = args + cmd2_args[arg] + str(funckeys[arg]) + ' '

    return args


def createGeoGrid(sub_path, input_geo_grid, run):
    for i in range(len(input_geo_grid)):
        fileloading.loadGeoModel(input_geo_grid[i], run)
        importlib.reload(params)

        new_lband = sub_path["GeoDoth5Path"] + input_geo_grid[i] + "Lensing" + ".h5"
        new_rtray = sub_path["GeoDoth5Path"] + input_geo_grid[i] + "RayTracing" + ".h5"
        print("Computation of {} Lensing Bands".format(input_geo_grid[i]))
        """Computation of the lensing bands______________________________________________________"""
        subprocess.run(['python3 ' + EZPaths.aartPath + '/lensingbands.py '], shell=True)

        spin_case = params.spin_case
        i_case = params.i_case

        fnbands = path + "LensingBands_a_%s_i_%s.h5" % (params.spin_case, params.i_case)

        print("Reading file: ", fnbands)

        h5f = h5py.File(fnbands, 'r')

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

        '''Analytical Ray-tracing______________________________________________________'''
        print("Analytical Raytracing of {}".format(input_geo_grid[i]))
        subprocess.run(['python3 ' + EZPaths.aartPath + '/raytracing.py '], shell=True)
        fnrays1 = path + "Rays_a_%s_i_%s.h5" % (spin_case, i_case)

        h5f = h5py.File(fnrays1, 'r')

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

        # Move lensing bands and raytracing bands
        subprocess.run(["mv " + fnbands + ' ' + new_lband], shell=True)
        subprocess.run(["mv " + fnrays1 + ' ' + new_rtray], shell=True)


line = "\n________________________\n"


def creatIntensityGrid(sub_path:dict, run:str, input_geo_grid_names:list[str], geo_grid_list,
                       intensity_models:list[tuple], var_params:list, constant_params:list,
                       total_models_count:int, action:dict):
    """

    Args:
        total_models_count:
        constant_params: parameters that will remain constant for all models
        var_params: list of key names to be varable parameters in the current run,
                    correspond to more than 1 entry in list for grid params
        sub_path: dictionary for all file paths
        input_geo_grid_names: list of the params files
        run: run name
        intensity_models: list[tuple(model_name,brightparams)]
        geo_grid_list: list[tuple(geo parameter, value)]
        action:

    Returns:

    """

    funckeys = {
        "emodelkey": 0,  # emodelkey Emission Model choice, 0 = thermal ultrarelativistic, 1 = power law
        "bkey": 2,  # bkey
        "nnoisykey": 0,  # nnoisykey Inoisy density. 0 = no noise, 1 = noise
        "tnoisykey": 0,  # tnoisykey Inoisy temperature
        "bnoisykey": 0  # bnoisykey Inoisy magnetic field
    }

    all_intent_names = []
    all_total_names = []
    all_bright_params = []
    all_230_total_jy_thin = []
    all_230_total_jy_thick = []

    print("All GeoGrid Names:  " + "\n" + str(input_geo_grid_names))
    for j in range(len(input_geo_grid_names)):
        for i in range(len(intensity_models)):
            current_geo_model = input_geo_grid_names[j]
            print("Current Geo Model = " + current_geo_model)
            fileloading.loadGeoModel(current_geo_model, run)
            lband = sub_path["GeoDoth5Path"] + current_geo_model + "Lensing" + ".h5"
            rtray = sub_path["GeoDoth5Path"] + current_geo_model + "RayTracing" + ".h5"

            # String Names
            all_intent_names += [intensity_models[i][0]]
            all_total_names += [current_geo_model + all_intent_names[i].replace("Model", "")]

            # ________________________________
            all_bright_params += [intensity_models[i][1]]
            # TODO: FIX
            print("Normalizing " + all_total_names[i])
            # all_bright_params[i]["n_th0"] = normalizingBrightparams.normalize(lband,rtray,all_bright_params[i])
            print("\n" + all_total_names[i] + "Normalized with a value of n_th0="
                  + str(all_bright_params[i]["n_th0"]) + "\n")

            print("Creating Intensity Movie for Model ", all_total_names[i])
            # intermodel_data = movieMakerIntensity.intensity_movie(
            #     action, sub_path, all_total_names[i], 2, all_bright_params[i])
            #
            # all_230_total_jy_thin += [intermodel_data["thin_total_flux"]]
            # all_230_total_jy_thick += [intermodel_data["thick_total_flux"]]


    # Make Docstring_____________________________________________

    doc_string_file = EZPaths.modelRunsDir + run + "/" + "AllModels.txt"
    cmd = "touch " + doc_string_file
    subprocess.run([cmd], shell=True)

    # Astrophysical part
    full_string = intensityModelsDocString(all_intent_names, all_bright_params,
                                           var_params, constant_params, total_models_count)

    # Geometrical Part
    geo_string = geoModelDocString(geo_grid_list, input_geo_grid_names)

    # writing
    doc_string_file = open(doc_string_file, 'w')
    # doc_string_file.write(full_string + geo_models_string)
    doc_string_file.write(full_string + geo_string)
    doc_string_file.close()

    brightparams_numpy_name = sub_path["meta"] + "AllBrightParamsList"
    all_full_names_numpy_name = sub_path["meta"] + "AllModelsList"
    all_230_total_jy_thin_numpy_name = sub_path["meta"] + "thin_total_flux"
    all_230_total_jy_thick_numpy_name = sub_path["meta"] + "thick_total_flux"

    np.save(all_full_names_numpy_name, np.array(all_total_names))
    np.save(brightparams_numpy_name, np.array(all_bright_params))
    np.save(all_230_total_jy_thin_numpy_name, np.array(all_230_total_jy_thin))
    np.save(all_230_total_jy_thick_numpy_name, np.array(all_230_total_jy_thick))


def geoModelDocString(geo_grid_list:list[tuple[list]], input_geo_grid_names):

    geo_models_string = "\nGEOMODELS\n"
    for i in range(len(geo_grid_list)):
        geo_models_string += input_geo_grid_names[i] + "| "
        for j in range(len(geo_grid_list[i][0])):
            geo_models_string += '\n    ' + geo_grid_list[i][0][j] + ": " + geo_grid_list[i][1][j]
        geo_models_string += '\n'

    return geo_models_string


def intensityModelsDocString(all_intent_names, all_bright_params, var_params, constant_params, total_models_count):
    """

    Args:
        all_intent_names: [intensity name]
        all_bright_params: [brightparams]
        var_params: list of key names to be varable parameters in the current run,
              correspond to more than 1 entry in list for grid params
        constant_params: list of parameters that will remain constant for all models
        total_models_count: total number of models for this run

    Returns:

    """

    line_small = "________________________________ \n"
    breaker = "     "

    string = line_small + line_small + line_small + "Total Number of Models: " + str(total_models_count) + '\n' + "Constant Params: " + '\n'
    for key in list(constant_params):
        if key != "n_th0":
            string += breaker + key + ": " + str(constant_params[key]) + '\n'

    string += line_small
    for i in range(len(all_intent_names)):
        current_name = all_intent_names[i]
        current_model = all_bright_params[i]
        string += line_small + current_name + '\n'

        for k in range(len(var_params)):
            string += breaker + var_params[k] + ": " + str(current_model[var_params[k]]) + '\n'
        string += breaker + "n_th0: " + str(current_model["n_th0"]) + '\n'

    return string + line_small + line_small + line_small


def graphCreation(sub_path, run, action, intent_grid_type=2):
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
    print(line)
    print(line)
    print("Initializing Graph Creation")
    all_full_model_names = np.load(sub_path["meta"] + "AllModelsList.npy")
    all_brightparams = np.load(sub_path["meta"] + "AllBrightParamsList.npy", allow_pickle=True)
    thin_total_flux = np.load(sub_path["meta"] + "thin_total_flux.npy")
    thick_total_flux = np.load(sub_path["meta"] + "thick_total_flux.npy")
    # current_intensity_model_name = intensity_models[i][0]
    # current_intensity_model = intensity_models[i][1]
    # full_current_model_name = current_geo_model + current_intensity_model_name.replace("Model", "")
    j = 0
    hist_flux_peaks_thins = []
    hist_flux_peaks_thicks = []
    hist_convs = []
    dim = [10, 8]
    for model in all_full_model_names:
        print(line)
        print("Running " + model)

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

        '''Reading of the lensing bands----------------------------------'''
        fnbands = lband

        print("Reading file: ", fnbands)

        h5f = h5py.File(fnbands, 'r')

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
        fnrays = rtray

        print("Reading file: ", fnrays)

        h5f = h5py.File(fnrays, 'r')

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

        num_of_intensity_points = janksys_thin[:,0].shape[0]
        print("Number of Intensity Points: ", num_of_intensity_points)

        xaxis = np.array(x_variable) / astroModels.scale_label[action['var']]
        # one_M = ilp.rg_func(brightparams["mass"] * u.g).to(u.m)
        # M2uas = np.arctan(one_M.value / dBH) / muas_to_rad

        fluxVNu_path = sub_path["fluxPath"] + model + "/"
        radVNu_path = sub_path["radPath"] + model + "/"
        image_path = sub_path["imagePath"] + model + "/"
        Optical_depth_path = sub_path["opticalDepth"] + model + "/"

        file_creation = [fluxVNu_path, radVNu_path,image_path,Optical_depth_path]

        for i in range(len(file_creation)):
            if not os.path.isdir(file_creation[i]):
                subprocess.run(["mkdir " + file_creation[i]], shell=True)

            else:
                subprocess.run(["rm -r " + file_creation[i]], shell=True)
                subprocess.run(["mkdir " + file_creation[i]], shell=True)

            print("Subdirectory '{}' created".format(file_creation[i]))

        # Points of Interest



        conv_1 = action["start"] + action["step"] * ilp.ring_convergance(mean_radii_Thick[:, 2], mean_radii_Thick[:, 3],
                                                                         3)
        conv_1 = conv_1 / astroModels.scale_label[action['var']]

        flux_peak_thin = action["start"] + action["step"] * np.argmax(janksys_thin[:, 3])
        flux_peak_thin = flux_peak_thin / astroModels.scale_label[action['var']]

        flux_peak_thick = action["start"] + action["step"] * np.argmax(janksys_thick[:, 3])
        flux_peak_thick = flux_peak_thick / astroModels.scale_label[action['var']]

        hist_flux_peaks_thins += [flux_peak_thin]
        hist_flux_peaks_thicks += [flux_peak_thick]
        hist_convs += [conv_1]

        poi = {
            "r_outer": r_outer,
            "flux_peak_thin": flux_peak_thin,
            "flux_peak_thick":flux_peak_thick,
            "conv_1": conv_1,
        }

        conv_1_style = {
            "color": 'dimgrey',
            "linestyle": "-",
            "linewidth": 2
        }

        r_outer_style = {
            "color": 'dimgrey',
            "linestyle": "-",
            "linewidth": 5
        }

        flux_peak_style = {
            "color": 'k',
            "linestyle": "-.",
            "linewidth": 3
        }

        # '''
        # "fluxPath"
        # "radPath"
        # "imagePath"
        # '''
        # _______________________________________________
        # _______________________________________________
        # ________________________________
        '''JANKSKY PLOTS----------------------------------'''

        fig, (ax, ax1) = plt.subplots(2, 1, figsize=dim, dpi=400, sharex=True)

        astroPloting.fluxThickThin(ax, ax1, xaxis, janksys_thin, janksys_thick,
                                   poi, conv_1_style, r_outer_style, flux_peak_style, action)

        figname = fluxVNu_path + model + "Flux.jpg"
        plt.savefig(figname, bbox_inches='tight')
        plt.close()
        print("Image '{}' Created".format(figname))

        # ______________________________________________
        # ______________________________________________
        # ___________________________________
        '''RADII PLOTS----------------------------------'''
        fig, (ax, ax1) = plt.subplots(2, 1, figsize=dim, dpi=400, sharex=True)

        astroPloting.radiiThickThin(ax, ax1, xaxis, mean_radii_Thin, mean_radii_Thick,
                                    poi, conv_1_style, r_outer_style, flux_peak_style, action)

        figname = radVNu_path + model + "Radii.jpeg"
        plt.savefig(figname, bbox_inches='tight')
        print("Image '{}' Created".format(figname))
        plt.close()

        # --------------------------------------------------
        # --------------------------------------------------
        # --------------------------------------------------
        '''Optical Depth----------------------------------'''
        fig, ax = plt.subplots(1, 1, figsize=dim, dpi=400)

        astroPloting.opticalDepth(ax, xaxis,
                                  [mean_optical_depth_I0,mean_optical_depth_I1,mean_optical_depth_I2],
                                  poi, conv_1_style, flux_peak_style, action)

        figname = Optical_depth_path + model + "OpticalDepth.jpeg"
        plt.savefig(figname, bbox_inches='tight')
        print("Image '{}' Created".format(figname))
        plt.close()

        '''Full Images----------------------------------'''
        k = action["start"]
        print("Constructing Full images for " + model)
        for i in range(num_of_intensity_points):
            brightparams = all_brightparams[j]
            brightparams["nu0"] = k
            print("Full image production for intensity frame: ", i)
            print(R"Observation frequency $\nu=$",k)
            # optically thin radii

            thin_alpha0 = radii_I0_Thin[i,:] * np.cos(theta)
            thin_beta0 = radii_I0_Thin[i,:] * np.sin(theta)
            thin_alpha1 = radii_I1_Thin[i,:] * np.cos(theta)
            thin_beta1 = radii_I1_Thin[i,:] * np.sin(theta)
            thin_alpha2 = radii_I2_Thin[i,:] * np.cos(theta)
            thin_beta2 = radii_I2_Thin[i,:] * np.sin(theta)
            thin_alpha_full = radii_Full_Thin[i,:] * np.cos(theta)
            thin_beta_full = radii_Full_Thin[i,:] * np.sin(theta)

            # full solution radii
            thick_alpha0 = radii_I0_Thick[i,:] * np.cos(theta)
            thick_beta0 = radii_I0_Thick[i,:] * np.sin(theta)
            thick_alpha1 = radii_I1_Thick[i,:] * np.cos(theta)
            thick_beta1 = radii_I1_Thick[i,:] * np.sin(theta)
            thick_alpha2 = radii_I2_Thick[i,:] * np.cos(theta)
            thick_beta2 = radii_I2_Thick[i,:] * np.sin(theta)
            thick_alpha_full = radii_FullAbsorption_Thick[i,:] * np.cos(theta)
            thick_beta_full = radii_FullAbsorption_Thick[i,:] * np.sin(theta)

            current_intensity_file = (sub_path["intensityPath"] + model + "/" + action["var"]
                                      + "_" + "{:.5e}".format(brightparams[action["var"]]))

            lim0 = 25

            print("Reading file: ", current_intensity_file)

            h5f = h5py.File(current_intensity_file, 'r')

            I0 = h5f['bghts0'][:]  # This implies I0 is 1 pass
            I1 = h5f['bghts1'][:]
            I2 = h5f['bghts2'][:]

            # I2_Absorb = h5f['bghts2_absorbtion'][:]
            # I1_Absorb = h5f['bghts1_absorbtion'][:]
            # I0_Absorb = h5f['bghts0_absorbtion'][:]
            Absorbtion_Image = h5f['bghts_full_absorbtion'][:]

            h5f.close()

            vmax0 = np.nanmax(I0 + I1 + I2) * 1.2
            fig, (ax0, ax1) = plt.subplots(1, 2, figsize=[15, 7], dpi=400)

            # Optically Thin

            im0 = ax0.imshow(I0 + I1 + I2, vmax=vmax0, origin="lower", cmap="afmhot", extent=[-lim0, lim0, -lim0, lim0])

            ax0.set_xlim(-10, 10)  # units of M
            ax0.set_ylim(-10, 10)

            ax0.text(-9, 8.5, astroModels.var_label[action["var"]]
                     + str(round(x_variable[i] / astroModels.scale_label[action["var"]], 2))
                     + ' ' + astroModels.units_label[action["var"]], fontsize=12, color="w")

            ax0.set_xlabel(r"$\alpha$" + " " + r"($\mu as$)")
            ax0.set_ylabel(r"$\beta$" + " " + r"($\mu as$)")
            ax0.title.set_text('Optically Thin Assumption')

            # Optically thick

            im1 = ax1.imshow(Absorbtion_Image, origin="lower", cmap="afmhot", extent=[-lim0, lim0, -lim0, lim0])

            #
            ax1.set_xlim(-10, 10)  # units of M
            ax1.set_ylim(-10, 10)

            ax1.set_xlabel(r"$\alpha$" + " " + r"($\mu as$)")

            ax1.title.set_text('Full Solution')

            colorbar0 = fig.colorbar(im1, fraction=0.046, pad=0.04, format=ticker.FuncFormatter(fmt), ticks=[
                vmax0 * .8,
                vmax0 * .6,
                vmax0 * .4,
                vmax0 * .2,
                vmax0 * .05
            ],
                                     ax=ax0
                                     )
            colorbar1 = fig.colorbar(im1, fraction=0.046, pad=0.04, format=ticker.FuncFormatter(fmt), ticks=[
                vmax0 * .8,
                vmax0 * .6,
                vmax0 * .4,
                vmax0 * .2,
                vmax0 * .05
            ],
                                     label="Brightnes Temperature (1e9 K)",
                                     ax=ax1
                                     )


            '''Radii Calc______________________'''
            # Thin
            lineCum_thickness = 4
            line0_thickness = 3
            line1_thickness = 2
            line2_thickness = 1

            ax0.plot(thin_alpha_full, thin_beta_full, color='tab:purple', linestyle='--', linewidth=lineCum_thickness)
            ax0.plot(thin_alpha0, thin_beta0, color='tab:red', linestyle='-', linewidth=line0_thickness)
            ax0.plot(thin_alpha1, thin_beta1, color='tab:orange', linestyle=':', linewidth=line1_thickness)
            ax0.plot(thin_alpha2, thin_beta2, color='tab:blue', linestyle='--', linewidth=line2_thickness)

            # Thick
            ax1.plot(thick_alpha_full, thick_beta_full, color='tab:purple', linestyle='--', linewidth=lineCum_thickness, label='Cumulative')
            ax1.plot(thick_alpha0, thick_beta0, color='tab:red', linestyle='-', linewidth=line0_thickness, label=R'n=0')
            ax1.plot(thick_alpha1, thick_beta1, color='tab:orange', linestyle=':', linewidth=line1_thickness, label=R'n=1')
            ax1.plot(thick_alpha2, thick_beta2, color='tab:blue', linestyle='--', linewidth=line2_thickness, label=R'n=2')
            ax1.legend()

            plt.subplots_adjust(wspace=.3)

            pltname = (image_path + 'FullImage_' + str(i) + "_Nu_"
                       + str(round(x_variable[i] / astroModels.scale_label[action["var"]], 2)) + ".jpeg")
            plt.savefig(pltname, bbox_inches='tight')
            print("Jpeg Created:  " + pltname)
            plt.close()

            # Get total jansky


            k += action['step']
        j += 1  # marker for which brightparams to use
    # histograms
    print(line)
    print("Creating Histograms")
    """Flux Peaks_____________________________________________________________________"""
    peak_hist_thin_path = sub_path["peakHistThin"]
    peak_hist_thick_path = sub_path["peakHistThick"]
    conv_hist_path = sub_path["convHist"]
    total_flux_path = sub_path["totalFlux"]

    # Thin_________________

    fig, ax = plt.subplots(1, 1, figsize=dim, dpi=400)

    astroPloting.histogram(ax,hist_flux_peaks_thins,"Flux Peak","Optically Thin Assumption Frequency")

    figname = peak_hist_thin_path + "FluxPeakThin.jpeg"
    plt.savefig(figname, bbox_inches='tight')
    print("Image '{}' Created".format(figname))
    plt.close()

    # Thick_______________

    fig, ax = plt.subplots(1, 1, figsize=dim, dpi=400)

    astroPloting.histogram(ax, hist_flux_peaks_thicks, "Flux Peak", "Full Solution Frequency")

    figname = peak_hist_thick_path + "FluxPeakThick.jpeg"
    plt.savefig(figname, bbox_inches='tight')
    print("Image '{}' Created".format(figname))
    plt.close()

    """Conv_____________________________________________________________________"""
    fig, ax = plt.subplots(1, 1, figsize=dim, dpi=400)

    astroPloting.histogram(ax, hist_convs, "Flux Peak", "Full Solution Frequency")

    figname = conv_hist_path + "convHist.jpeg"
    plt.savefig(figname, bbox_inches='tight')
    print("Image '{}' Created".format(figname))
    plt.close()

    """230 total flux Thin___________________________________________________________"""
    fig, ax = plt.subplots(1, 1, figsize=dim, dpi=400)

    astroPloting.bar(ax, np.arange(len(all_full_model_names)),thin_total_flux, "Total Flux at 230GHz",
                     "Optically Thin Assumption Frequency",all_full_model_names)

    figname = total_flux_path + "thin.jpeg"
    plt.savefig(figname, bbox_inches='tight')
    print("Image '{}' Created".format(figname))
    plt.close()

    """230 total flux Thick___________________________________________________________"""
    fig, ax = plt.subplots(1, 1, figsize=dim, dpi=400)

    astroPloting.bar(ax, np.arange(len(all_full_model_names)), thick_total_flux, "Total Flux at 230GHz",
                     "Full Solution Frequency", all_full_model_names)

    figname = total_flux_path + "thick.jpeg"
    plt.savefig(figname, bbox_inches='tight')
    print("Image '{}' Created".format(figname))
    plt.close()

def fmt(x, pos):
    x = x / 1e9
    return '{:.2f}'.format(x)

# Full Run
# current_run = "run1"
# current_bp = astroModels.bp_run1
# current_var_params = ["p_temp", "p_mag"]
# current_geo_grid = ["ModelA", "ModelB"]
