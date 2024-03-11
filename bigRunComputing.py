import subprocess
import EZPaths
import os
from aart_func import *
from params import *
import importlib
import params
import astroModels
import fileloading
from movieMakerV2 import movieMakerIntensity

action = {
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


def creatIntensityGrid(sub_path,input_geo_grid,run,intensity_models,full_string,geo_grid_values):
    funckeys = {
        "emodelkey": 0,  # emodelkey Emission Model choice, 0 = thermal ultrarelativistic, 1 = power law
        "bkey": 2,  # bkey
        "nnoisykey": 0,  # nnoisykey Inoisy density. 0 = no noise, 1 = noise
        "tnoisykey": 0,  # tnoisykey Inoisy temperature
        "bnoisykey": 0  # bnoisykey Inoisy magnetic field
    }

    all_intent_names = []

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
                                                full_current_model_name,2,current_intensity_model)

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
            # all_intent_names += [full_current_model_name]

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
    numpy_name = EZPaths.modelRunsDir + run + "/" + "AllModelsList"
    np.save(numpy_name, all_intent_names)


def creatDataGrids(sub_path,run,intensity_models):
    """
        sub_paths = {
        "GeoDoth5Path"
        "intensityPath"
        "fluxPath"
        "radPath"
        "imagePath"

    """
    run_path = EZPaths.modelRunsDir + run + "/"
    all_intent_names = np.load(run_path + "AllModelsList.npy")

    for model in all_intent_names:
        intensity_file = sub_path["intensityPath"] + model + ".h5"
        print("Reading file: ", intensity_file)

        h5f = h5py.File(intensity_file, 'r')

        I0 = h5f['bghts0'][:]  # This implies I0 is 1 pass
        I1 = h5f['bghts1'][:]
        I2 = h5f['bghts2'][:]

        I2_Absorb = h5f['bghts2_absorbtion'][:]
        I1_Absorb = h5f['bghts1_absorbtion'][:]
        I0_Absorb = h5f['bghts0_absorbtion'][:]
        Absorbtion_Image = h5f['bghts_full_absorbtion'][:]

        h5f.close()

        thin_assumption = [I0, I1, I2, I0 + I1 + I2]
        full_solution = [I0_Absorb, I1_Absorb, I2_Absorb, Absorbtion_Image]
        rad_v_flux_file = sub_path["fluxPath"] + model + "/"
        isDir = os.path.exists(rad_v_flux_file)
        if not isDir:
            os.makedirs(rad_v_flux_file)
            print("Flux Versus Subdirectory {} Created".format(rad_v_flux_file))



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




