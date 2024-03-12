import subprocess

from matplotlib import ticker

import EZPaths
import os
from aart_func import *
from params import *
import importlib
import params
import astroModels
import fileloading
from movieMakerV2 import movieMakerIntensity
from astropy import units as u

kw_action = {
    "var": "nu0",
    "start": 3.00e+10,
    "stop": 6.00e+10,
    "step": 1.00e+10
}


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

    funckeys = {
        "emodelkey": 0,  # emodelkey Emission Model choice, 0 = thermal ultrarelativistic, 1 = power law
        "bkey": 2,  # bkey
        "nnoisykey": 0,  # nnoisykey Inoisy density. 0 = no noise, 1 = noise
        "tnoisykey": 0,  # tnoisykey Inoisy temperature
        "bnoisykey": 0  # bnoisykey Inoisy magnetic field
    }

    # brightparams = fpp.bp_steeperT
    # funckeys = fpp.fk_fiducial

    for arg in cmd1_args:
        args = args + cmd1_args[arg] + str(brightparams[arg]) + ' '

    for arg in cmd2_args:
        args = args + cmd2_args[arg] + str(funckeys[arg]) + ' '

    return args


def createGeoGrid(sub_path,input_geo_grid,run):
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


def creatIntensityGrid(sub_path,input_geo_grid,run,intensity_models,full_string,geo_grid_values,action):
    funckeys = {
        "emodelkey": 0,  # emodelkey Emission Model choice, 0 = thermal ultrarelativistic, 1 = power law
        "bkey": 2,  # bkey
        "nnoisykey": 0,  # nnoisykey Inoisy density. 0 = no noise, 1 = noise
        "tnoisykey": 0,  # tnoisykey Inoisy temperature
        "bnoisykey": 0  # bnoisykey Inoisy magnetic field
    }

    all_intent_names = []
    all_bright_params = []
    for j in range(len(input_geo_grid)):
        for i in range(len(intensity_models)):
            current_geo_model = input_geo_grid[j]
            fileloading.loadGeoModel(current_geo_model,run)
            lband = sub_path["GeoDoth5Path"] + current_geo_model + "Lensing" + ".h5"
            rtray = sub_path["GeoDoth5Path"] + current_geo_model + "RayTracing" + ".h5"

            # Intensity

            current_intensity_model_name = intensity_models[i][0]
            current_intensity_model = intensity_models[i][1]
            full_current_model_name = current_geo_model + current_intensity_model_name.replace("Model","")

            movieMakerIntensity.intensity_movie(action, sub_path,
                                                full_current_model_name, 2, current_intensity_model)

            all_intent_names += [full_current_model_name]
            all_bright_params += [current_intensity_model]
            # args = createIntensityArgs(current_intensity_model)
            #
            # args += "--lband " + lband + " --rtray " + rtray
            #
            # subprocess.run(['python3 ' + EZPaths.aartPath + '/radialintensity.py' + args], shell=True)
            #
            # new_intensity_path = sub_path["intensityPath"] + full_current_model_name + "Intensity" + ".h5"
            #
            # fnrays = fileloading.intensityNameNoUnits(current_intensity_model, funckeys)
            #
            # subprocess.run(["mv " + fnrays + ' ' + new_intensity_path], shell=True)
            #

    # Make Docstring
    doc_string_file = EZPaths.modelRunsDir + run + "/" + "AllModels.txt"
    cmd = "touch " + doc_string_file
    subprocess.run([cmd], shell=True)

    geo_models_string = "\nGEOMODELS\n"
    for i in range(len(geo_grid_values)):
        geo_models_string += input_geo_grid[i] + "| "
        for j in range(len(geo_grid_values[i][0])):
            geo_models_string += '\n    ' + geo_grid_values[i][0][j] + ": " + geo_grid_values[i][1][j]
        geo_models_string += '\n'

    doc_string_file = open(doc_string_file,'w')
    doc_string_file.write(full_string + geo_models_string)
    doc_string_file.close()

    all_intent_names=np.array(all_intent_names)
    all_bright_params = np.array(all_bright_params)

    bright_numpy_name = EZPaths.modelRunsDir + run + "/" + "AllBrightParamsList"
    numpy_name = EZPaths.modelRunsDir + run + "/" + "AllModelsList"

    np.save(numpy_name, all_intent_names)
    np.save(bright_numpy_name,all_bright_params)


def graphCreation(sub_path,run,action,intent_grid_type=2):
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

    all_intent_names = np.load(EZPaths.modelRunsDir + run + "/" + "AllModelsList.npy")
    all_brightparams = np.load(EZPaths.modelRunsDir + run + "/" + "AllBrightParamsList.npy", allow_pickle=True)
    # current_intensity_model_name = intensity_models[i][0]
    # current_intensity_model = intensity_models[i][1]
    # full_current_model_name = current_geo_model + current_intensity_model_name.replace("Model", "")

    scale_label = {
        'nu0': 1e9,  # GHz
        'mass': 1e9 * 1.989e33,  # Billion Solar Masses
        'scale_height': 1,  # Rg
        'theta_b': 1,  # Rads
        'beta': 1,
        'r_ie': 1,
        'rb_0': 1,  # Rg
        'n_th0': 1 / 1e6,  # 1/cm^3
        't_e0': 1e9,  # GK
        'p_dens': 1,
        'p_temp': 1,
        'nscale': 1
    }

    for model in all_intent_names:
        print(line)
        print(line)
        print("Running " + model)
        brightparams = all_brightparams

        current_geo_model = model[0:len(model) - intent_grid_type]
        fileloading.loadGeoModel(current_geo_model, run)
        lband = sub_path["GeoDoth5Path"] + current_geo_model + "Lensing" + ".h5"
        rtray = sub_path["GeoDoth5Path"] + current_geo_model + "RayTracing" + ".h5"

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
        radii_FullAbsorption_Thick = np.load(data_path + "radii_FullAbsorption_Thick.npy")
        radii_I0_Thick = np.load(data_path + "radii_I0_Thick.npy")
        radii_I1_Thick = np.load(data_path + "radii_I1_Thick.npy")
        radii_I2_Thick = np.load(data_path + "radii_I2_Thick.npy")
        theta = np.load(data_path + "theta.npy")

        dim = [10, 8]
        xaxis = np.array(x_variable) / scale_label[action['var']]
        # one_M = ilp.rg_func(brightparams["mass"] * u.g).to(u.m)
        # M2uas = np.arctan(one_M.value / dBH) / muas_to_rad

        fluxVNu_path = sub_path["fluxPath"] + model + "/"
        radVNu_path = sub_path["radPath"] + model + "/"
        image_path = sub_path["imagePath"] + model + "/"

        if not os.path.isdir(fluxVNu_path):
            subprocess.run(["mkdir " + fluxVNu_path], shell=True)



        # '''
        # "fluxPath"
        # "radPath"
        # "imagePath"
        # '''

        # JANKSKY PLOTS________________________________

        fig, (ax,ax1) = plt.subplots(2, 1, figsize=dim, dpi=400)
        ax.plot(xaxis, janksys_thin[:, 0], '-', label='n=0', color='tab:red', linewidth=3)
        ax.plot(xaxis, janksys_thin[:, 1], ':', label='n=1', color='tab:orange', linewidth=3)
        ax.plot(xaxis, janksys_thin[:, 2], '--', label='n=2', color='tab:blue', linewidth=3)
        ax.plot(xaxis, janksys_thin[:, 3], '-.', label='Total', color='tab:purple', linewidth=3)

        flux_peak = action["start"] + action["step"] * np.argmax(janksys_thin[:, 3])
        flux_peak = flux_peak / scale_label[action['var']]

        ax.axhline(.5, color='k', label=R'.5 $J_y$', linestyle=":")
        ax.axvline(230, color='k', linestyle=":")

        # Labels
        ax.set_ylabel("Total Flux ({})".format(R'$J_y$'))
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.xaxis.set_minor_formatter(ticker.FormatStrFormatter('%.1f'))
        ax.xaxis.set_major_formatter(ticker.FormatStrFormatter("%.1f"))
        ax.yaxis.set_minor_formatter(ticker.FormatStrFormatter('%.1f'))
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.1f"))

        ax.tick_params('x', which="both", labelbottom=False)
        ax.tick_params('y', which="minor", labelleft=False)
        n = 4  # Keeps every 4th label
        [l.set_visible(False) for (i, l) in enumerate(ax.xaxis.get_minorticklabels()) if i % n != 0]
        ax.tick_params('both', length=10, width=1, which='major')
        ax.set_xlim(xaxis[0], xaxis[xaxis.size - 1])
        ax.legend(loc='lower left')

        conv_1 = action["start"] + action["step"] * ilp.ring_convergance(mean_radii_Thick[:, 2], mean_radii_Thick[:, 3],
                                                                         5)
        conv_1 = conv_1 / scale_label[action['var']]

        ax1.plot(xaxis, janksys_thick[:, 0], '-', label='One pass', color='tab:red', linewidth=3)
        ax1.plot(xaxis, janksys_thick[:, 1], ':', label='Two passes', color='tab:orange', linewidth=3)
        ax1.plot(xaxis, janksys_thick[:, 2], '--', label='Three passes', color='tab:blue', linewidth=3)
        ax1.plot(xaxis, janksys_thick[:, 3], '-.', label='Cumulative', color='tab:purple', linewidth=3)

        flux_peak = action["start"] + action["step"] * np.argmax(janksys_thick[:, 3])
        flux_peak = flux_peak / scale_label[action['var']]

        ax1.axhline(.5, color='k', label=R'.5 $J_y$', linestyle=":")
        ax1.axvline(230, color='k', linestyle=":")
        ax1.axvline(conv_1, color='k', linestyle="--", linewidth=3)
        ax1.axvline(flux_peak, color='dimgrey', linestyle="-", linewidth=2)

        # Labels
        ax1.set_ylabel("Total Flux ({})".format(R'$J_y$'))
        ax1.set_xscale('log')
        ax1.set_yscale('log')

        ax1.xaxis.set_minor_formatter(ticker.FormatStrFormatter('%.1f'))
        ax1.xaxis.set_major_formatter(ticker.FormatStrFormatter("%.1f"))
        # ax1.yaxis.set_minor_formatter(ticker.FormatStrFormatter('%.4f'))
        ax1.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.0e"))

        ax1.tick_params('x', which="both", labelbottom=False)
        ax1.tick_params('y', which="minor", labelleft=False)

        n = 4  # Keeps every 4th label
        [l.set_visible(False) for (i, l) in enumerate(ax1.xaxis.get_minorticklabels()) if i % n != 0]
        ax1.tick_params('both', length=10, width=1, which='major')
        ax1.set_xlim(xaxis[0], xaxis[xaxis.size - 1])

        ax1.legend(loc='lower left')

        plt.savefig(fluxVNu_path + "flux.jpg")
        plt.close()
        #
        # # RADII PLOTS___________________________________
        # ax[1] = plt.subplot(2, 1, 2, sharex=ax[0])
        # # ax[1].axhline(r_inner, color='k', linewidth=3, linestyle=":")  # , label='Blackhole Inner Shadow'
        # ax[1].axhline(r_outer, color='dimgrey', linewidth=5)  # , label='Blackhole Outer Shadow'
        # ax[1].plot(xaxis, mean_radii_Thin[:, 0], '-', label='n=0', color='tab:red', linewidth=3)
        # ax[1].plot(xaxis, mean_radii_Thin[:, 1], ':', label='n=1', color='tab:orange', linewidth=3)
        # ax[1].plot(xaxis, mean_radii_Thin[:, 2], '-.', label='n=2', color='tab:blue', linewidth=3)
        #
        # ax[1].axvline(flux_peak, color='k', linestyle="-.")
        # ax[1].axvline(230, color='k', linestyle=":")
        # # Labels
        #
        # ax[1].set_xscale('log')
        # ax[1].set_yscale('log')
        #
        # ax[1].set_xlabel(var_label[action["var"]].replace('=', '') + ' (' + units_label[action["var"]] + ')')
        # ax[1].set_ylabel("Ring Radii ({})".format(R'$R_g$'))
        #
        # ax[1].xaxis.set_minor_formatter(ticker.FormatStrFormatter('%.0f'))
        # ax[1].xaxis.set_major_formatter(ticker.FormatStrFormatter("%.0f"))
        # ax[1].yaxis.set_minor_formatter(ticker.FormatStrFormatter('%.1f'))
        # ax[1].yaxis.set_major_formatter(ticker.FormatStrFormatter("%.1f"))
        #
        # n = 4  # Keeps every 4th label
        # [l.set_visible(False) for (i, l) in enumerate(ax[1].xaxis.get_minorticklabels()) if i % n != 0]
        #
        # ax[1].legend(frameon=False)
        # ax[1].set_xlim(xaxis[0], xaxis[xaxis.size - 1])
        #
        # new_ticks = [xaxis[0], 230, xaxis[xaxis.size - 1]]
        # ax[1].set_xticks(new_ticks)
        # ax[1].tick_params('x', length=20, width=1, which='major', labelrotation=90)
        # plt.savefig(final_graph_path + "Thin_" + iteration + ".png")
        # # Markers
        # plt.close()
        #
        # # Thick Full Image--------------------------------------------------------------------------------------------------
        #
        # conv_1 = action["start"] + action["step"] * ilp.ring_convergance(mean_radii_Thick[:, 2], mean_radii_Thick[:, 3],
        #                                                                  5)
        # conv_1 = conv_1 / scale_label[action['var']]
        #
        # fig = plt.subplots(2, 1, figsize=dim, dpi=400)
        # ax1 = [None, None]
        # ax1[0] = plt.subplot(2, 1, 1)
        # ax1[0].plot(xaxis, janksys_thick[:, 0], '-', label='One pass', color='tab:red', linewidth=3)
        # ax1[0].plot(xaxis, janksys_thick[:, 1], ':', label='Two passes', color='tab:orange', linewidth=3)
        # ax1[0].plot(xaxis, janksys_thick[:, 2], '--', label='Three passes', color='tab:blue', linewidth=3)
        # ax1[0].plot(xaxis, janksys_thick[:, 3], '-.', label='Cumulative', color='tab:purple', linewidth=3)
        #
        # flux_peak = action["start"] + action["step"] * np.argmax(janksys_thick[:, 3])
        # flux_peak = flux_peak / scale_label[action['var']]
        #
        # ax1[0].axhline(.5, color='k', label=R'.5 $J_y$', linestyle=":")
        # ax1[0].axvline(230, color='k', linestyle=":")
        # ax1[0].axvline(conv_1, color='k', linestyle="--", linewidth=3)
        # ax1[0].axvline(flux_peak, color='dimgrey', linestyle="-", linewidth=2)
        #
        # # Labels
        # ax1[0].set_ylabel("Total Flux ({})".format(R'$J_y$'))
        # ax1[0].set_xscale('log')
        # ax1[0].set_yscale('log')
        #
        # ax1[0].xaxis.set_minor_formatter(ticker.FormatStrFormatter('%.1f'))
        # ax1[0].xaxis.set_major_formatter(ticker.FormatStrFormatter("%.1f"))
        # # ax1[0].yaxis.set_minor_formatter(ticker.FormatStrFormatter('%.4f'))
        # ax1[0].yaxis.set_major_formatter(ticker.FormatStrFormatter("%.0e"))
        #
        # ax1[0].tick_params('x', which="both", labelbottom=False)
        # ax1[0].tick_params('y', which="minor", labelleft=False)
        #
        # n = 4  # Keeps every 4th label
        # [l.set_visible(False) for (i, l) in enumerate(ax1[0].xaxis.get_minorticklabels()) if i % n != 0]
        # ax1[0].tick_params('both', length=10, width=1, which='major')
        # ax1[0].set_xlim(xaxis[0], xaxis[xaxis.size - 1])
        #
        # ax1[0].legend(loc='lower left')
        #
        # # RADII PLOTS___________________________________
        # ax1[1] = plt.subplot(2, 1, 2, sharex=ax1[0])
        #
        # # ax1[1].axhline(r_inner, color='k', linewidth=2, linestyle=":")  # , label='Blackhole Inner Shadow'
        # ax1[1].axhline(r_outer, color='dimgrey', linewidth=5)  # , label='Blackhole Outer Shadow'
        # # ax1[1].scatter(xaxis[0],r_outer, marker="o", linewidth=10)
        # # ax1[1].scatter(xaxis[0],r_inner, marker="o", linewidth=10)
        #
        # ax1[1].plot(xaxis, mean_radii_Thick[:, 0], '-', label='One pass', color='tab:red', linewidth=3)
        # ax1[1].plot(xaxis, mean_radii_Thick[:, 1], ':', label='Two passes', color='tab:orange', linewidth=3)
        # ax1[1].plot(xaxis, mean_radii_Thick[:, 2], '--', label='Three passes', color='tab:blue', linewidth=3)
        # ax1[1].plot(xaxis, mean_radii_Thick[:, 3], '-.', label='Cumulative', color='tab:purple', linewidth=3)
        #
        # ax1[1].axvline(230, color='k', linestyle=":")
        # ax1[1].axvline(conv_1, color='k', linestyle="--", linewidth=3)
        # ax1[1].axvline(flux_peak, color='dimgrey', linestyle="-", linewidth=2)
        #
        # # ax1[1].axvline(conv_1, color='b',label=R'1% diff', linestyle="--")
        # # ax1[1].axvline(conv_5, color='g',label=R'5% diff', linestyle=":")
        # # ax1[1].axvline(conv_10, color='r',label=R'10% diff', linestyle="-.")
        #
        # # Labels
        # ax1[1].set_xlabel(var_label[action["var"]].replace('=', '') + ' (' + units_label[action["var"]] + ')')
        # ax1[1].set_ylabel("Ring Radii ({})".format(R'$R_g$'))
        # ax1[1].set_xscale('log')
        # ax1[1].set_yscale('log')
        #
        # ax1[1].xaxis.set_minor_formatter(ticker.FormatStrFormatter('%.0f'))
        # ax1[1].xaxis.set_major_formatter(ticker.FormatStrFormatter("%.0f"))
        # ax1[1].yaxis.set_minor_formatter(ticker.FormatStrFormatter('%.0f'))
        # ax1[1].yaxis.set_major_formatter(ticker.FormatStrFormatter("%.0f"))
        #
        # new_ticks = [xaxis[0], 230, conv_1, flux_peak, xaxis[xaxis.size - 1]]
        # ax1[1].set_xticks(new_ticks)
        #
        # # new_ticks = np.append(ax1[1].get_yticks(), r_outer)
        # # new_ticks = np.append(new_ticks, r_inner)
        # # print(new_ticks)
        # # ax1[1].set_yticks(new_ticks)
        # n = 4  # Keeps every 4th label
        # [l.set_visible(False) for (i, l) in enumerate(ax1[1].xaxis.get_minorticklabels()) if i % n != 0]
        # ax1[1].legend(frameon=False)
        # ax1[1].set_xlim(xaxis[0], xaxis[xaxis.size - 1])
        #
        # ax1[1].tick_params('x', length=20, width=1, which='major', labelrotation=80)
        #
        # plt.savefig(final_graph_path + "Thick_" + iteration + ".png")
        # plt.close()

#
# dim = [14,10]


# def createRadVFlux(rad_v_flux_file,thin_assumption,full_solution):
#     fig, ax = plt.subplots(1, 1, figsize=dim, dpi=400)
#     ax = plt.subplot(2, 1, 1)
#     ax.plot(xaxis, janksys_thin[:, 0], '-', label='n=0', color='tab:red', linewidth=3)
#     ax.plot(xaxis, janksys_thin[:, 1], ':', label='n=1', color='tab:orange', linewidth=3)
#     ax.plot(xaxis, janksys_thin[:, 2], '--', label='n=2', color='tab:blue', linewidth=3)
#     ax.plot(xaxis, janksys_thin[:, 3], '-.', label='Total', color='tab:purple', linewidth=3)
#
#     flux_peak = action["start"] + action["step"] * np.argmax(janksys_thin[:, 3])
#     flux_peak = flux_peak / scale_label[action['var']]
#
#     ax[0].axhline(.5, color='k', label=R'.5 $J_y$', linestyle=":")
#     ax[0].axvline(230, color='k', linestyle=":")
#
#     # Labels
#     ax[0].set_ylabel("Total Flux ({})".format(R'$J_y$'))
#     ax[0].set_xscale('log')
#     ax[0].set_yscale('log')
#     ax[0].xaxis.set_minor_formatter(ticker.FormatStrFormatter('%.1f'))
#     ax[0].xaxis.set_major_formatter(ticker.FormatStrFormatter("%.1f"))
#     ax[0].yaxis.set_minor_formatter(ticker.FormatStrFormatter('%.1f'))
#     ax[0].yaxis.set_major_formatter(ticker.FormatStrFormatter("%.1f"))
#
#     ax[0].tick_params('x', which="both", labelbottom=False)
#     ax[0].tick_params('y', which="minor", labelleft=False)
#     n = 4  # Keeps every 4th label
#     [l.set_visible(False) for (i, l) in enumerate(ax[0].xaxis.get_minorticklabels()) if i % n != 0]
#     ax[0].tick_params('both', length=10, width=1, which='major')
#     ax[0].set_xlim(xaxis[0], xaxis[xaxis.size - 1])
#     ax[0].legend(loc='lower left')


# Full Run
# current_run = "run1"
# current_bp = astroModels.bp_run1
# current_var_params = ["p_temp", "p_mag"]
# current_geo_grid = ["ModelA", "ModelB"]




