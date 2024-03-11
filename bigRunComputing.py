import subprocess
import EZPaths
import os
from aart_func import *
from params import *
import astroModels
import fileloading


def createGeoGrid(input_geo_grid, run, sub_path):
    for i in range(len(input_geo_grid)):
        fileloading.loadGeoModel(input_geo_grid[i], run)

        new_lband = sub_path["GeoDoth5Path"] + input_geo_grid[i] + "Lensing" + ".h5"
        new_rtray = sub_path["GeoDoth5Path"] + input_geo_grid[i] + "RayTracing" + ".h5"
        print("Computation of {} Lensing Bands".format(input_geo_grid[i]))
        """Computation of the lensing bands______________________________________________________"""
        subprocess.run(['python3 ' + EZPaths.aartPath + '/lensingbands.py '], shell=True)

        fnbands = path + "LensingBands_a_%s_i_%s.h5" % (spin_case, i_case)

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


def creatIntensityGrid(sub_path,input_geo_grid, intensity_models, run):
    funckeys = {
        "emodelkey": 0,  # emodelkey Emission Model choice, 0 = thermal ultrarelativistic, 1 = power law
        "bkey": 2,  # bkey
        "nnoisykey": 0,  # nnoisykey Inoisy density. 0 = no noise, 1 = noise
        "tnoisykey": 0,  # tnoisykey Inoisy temperature
        "bnoisykey": 0  # bnoisykey Inoisy magnetic field
    }

    for j in range(len(input_geo_grid)):
        for i in range(len(intensity_models)):
            current_geo_model = input_geo_grid[j]
            fileloading.loadGeoModel(current_geo_model,run)
            lband = sub_path["GeoDoth5Path"] + current_geo_model + "Lensing" + ".h5"
            rtray = sub_path["GeoDoth5Path"] + current_geo_model + "RayTracing" + ".h5"

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

            # Intensity

            current_intensity_model_name = intensity_models[i][0]
            current_intensity_model = intensity_models[i][1]

            full_current_model_name = current_geo_model + current_intensity_model_name.replace("Model","")

            args = createIntensityArgs(current_intensity_model)

            subprocess.run(['python3 ' + EZPaths.aartPath + '/radialintensity.py' + args], shell=True)

            new_intensity_path = sub_path["intensityPath"] + full_current_model_name + "RayTracing" + ".h5"

            fnrays = fileloading.intensityNameWrite(current_intensity_model, funckeys)

            subprocess.run(["mv " + fnrays + ' ' + new_intensity_path], shell=True)


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


geo_grid = ["ModelA", "ModelB"]
sub_paths, all_models, model_name_string = fileloading.runsInit("run1",astroModels.bp_run1, ["p_temp","p_mag"])

