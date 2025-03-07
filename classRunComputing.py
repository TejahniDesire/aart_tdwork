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
import astroModels
from scipy.interpolate import interp1d
import importlib
import params
from astroModels import *
import fileloading
from movieMakerV2 import movieMakerIntensity
import normalizingBrightparams
from astropy import units as u
from functools import partial 
import re

breaker = "     "
line = "\n________________________\n"
long_line = "\n________________________________________________________________________\n"
line_small = "________________________________ \n"
"/scratch/gpfs/td6241/aart/rawResults/Intensity_a_0.001_i_17_nu_2.30000e+11_mass_1.32497e+43_scaleh_0.5_thetab_0.873_beta_1.00_rie_10.0_rb_5.0_nth0_1.3e+05_te0_5.0e+10_b0_8.000e+00_pdens_-0.7_ptemp_-1.5_pmag_-2.0_nscale_0.4_emkey_0_bkey_2_nkey_0_tnkey_0_bnkey_0_magkey_0.h5"
"/scratch/gpfs/td6241/aart/rawResults/Intensity_a_0.9375_i_17_nu_2.30000e+11_mass_1.32497e+43_scaleh_0.5_thetab_0.873_beta_1.00_rie_10.0_rb_2.0_nth0_1.9e+04_te0_2.0e+11_b0_8.131e+00_pdens_-0.7_ptemp_-0.84_pmag_-0.3_nscale_0.4_emkey_0_bkey_2_nkey_0_tnkey_0_bnkey_0_magkey_0.h5"


class BigRuns:

    def __init__(self, run: str, intensity_grid_params: dict, var_intensity_grid_names, geo_grid_params, geo_grid_names,
                 normalized_brightparams=False,funcKey=None):
        """

        Args:
            run:
            intensity_grid_params: dictionary containing all intensity parameters. Each value is a
                         list containing all parameters to be computed in current run
            var_intensity_grid_names: list of key names to be varable parameters in the current run,
                        correspond to more than 1 entry in list for grid params
            geo_grid_params: list[tuple(geo parameter, value)]
            geo_grid_names: list of key names to be varable parameters in the current run,
                        correspond to more than 1 entry in list for grid params
        """

        """
        Args:
            run: name of this run
            intensity_grid_params: dictionary containing all intensity parameters. Each value is a
                         list containing all parameters to be computed in current run
            var_intensity_grid_names: list of key names to be varable parameters in the current run,
                        correspond to more than 1 entry in list for grid params
            geo_grid_list: list[tuple(geo parameter, value)]
        Returns:

        """
        self.intensity_grid_params = intensity_grid_params
        self.run = run
        if funcKey is None:
            self.funckeys = funckeys
        else:
            self.funckeys = funcKey
        self.var_intensity_grid_names = var_intensity_grid_names
        self.geo_grid_names = geo_grid_names
        self.geo_grid_params = geo_grid_params
        self.already_normalized_brightparams = normalized_brightparams

        # Create Directories______________________________________________________________________
        fileloading.creatSubDirectory(EZPaths.modelRunsDir, "storing all run results")
        main_path = EZPaths.modelRunsDir + run + "/"
        fileloading.creatSubDirectory(main_path, "storing this run result")

        image_path = main_path + "Images/"
        data_path = main_path + "Data/"
        meta_path = main_path + "InterModel/"
        fileloading.creatSubDirectory(image_path, "images")
        fileloading.creatSubDirectory(data_path, "data")
        fileloading.creatSubDirectory(meta_path, "inter model data")

        self.sub_paths = {
            "GeoDoth5Path": data_path + "geo/",
            "intensityPath": data_path + "intensity/",
            "runWideNumpy": data_path + "numpy/",
            "fluxPath": image_path + "fluxVNu/",
            "radPath": image_path + "radVNu/",
            "imagePath": image_path + "image/",
            "opticalDepth": image_path + "opticalDepth/",
            "peakHistThin": image_path + "peakHistThin/",
            "peakHistThick": image_path + "peakHistThick/",
            "convHist": image_path + "convHist/",
            "totalFlux": image_path + "totalFlux230Ghz/",
            "RadVVarphi": image_path + "radVVarphi/",
            "fluxVRadii": image_path + "fluxVVarphi/",
            "meta": meta_path,
            "3d": image_path + "3dPloting/",
            "movie": image_path + "Movie/"

        }

        for key in list(self.sub_paths):
            fileloading.creatSubDirectory(self.sub_paths[key])

        # Analyze run type______________________________________________________________________________________________
        self.total_intensity_models_count = 1  # total number of intensity models for this run
        self.run_type = 0  # number of varied intensity paramters
        self.variable_param_ranges = {}  # number of individual points in [key] paramter dimension
        self.constant_params = {}  # dict[key] = value of parameter [key] to be held constant

        for key in list(self.intensity_grid_params):

            param_range = len(self.intensity_grid_params[key])
            # print("Key, " + key + " count: " + str(var_params.count(key)))
            if param_range > 1 and self.var_intensity_grid_names.count(key) == 0:
                raise ValueError("More than one parameter put in for static parameter")

            if self.var_intensity_grid_names.count(key) == 1:
                if param_range <= 1:
                    raise ValueError("Static parameter listed as variable parameter")

                self.run_type += 1
                self.variable_param_ranges[key] = param_range
                self.total_intensity_models_count = self.total_intensity_models_count * param_range
            else:
                self.constant_params[key] = self.intensity_grid_params[key][0]

        print("Final Grid Type: ", self.run_type)
        # Create astro ParamGrid________________________________________________________________________________________
        self.string_order = {
            "p_mag":"{:.2f}",
            "p_temp":"{:.2f}",
            "p_dens":"{:.2f}",
            "theta_b":"{:.4f}",
            "mass":"{:.4e}",
            "t_e0":"{:.2e}",
            "nu0": "{:.2}",
            "scale_height":"{:.2f}",
            "rb_0":"{:.2f}",
            "beta":"{:.2f}",
            "b_0":"{:.2f}",
            "r_ie":"{:.2f}",
            "nscale":"{:.2f}",
        }
        if self.already_normalized_brightparams:

            file_paths = [
                self.sub_paths["runWideNumpy"] + "all_intensity_model_brightparams.npy",
                self.sub_paths["runWideNumpy"] + "all_intensity_model_names.npy",
                self.sub_paths["runWideNumpy"] + "all_model_brightparams.npy",
                self.sub_paths["runWideNumpy"] + "all_model_names.npy",
                self.sub_paths["runWideNumpy"] + "total_models_count.npy"
            ]
            all_intensity_model_brightparams = []
            all_intensity_model_names = []
            all_model_brightparams = []
            all_model_names = []
            total_models_count = 0

            arrays = [
                all_intensity_model_brightparams,
                all_intensity_model_names,
                all_model_brightparams,
                all_model_names,
                total_models_count
            ]
            k = 0
            for file in file_paths:
                arrays[k] = np.load(file, allow_pickle=True)
                k += 1

            self.all_intensity_model_brightparams = arrays[0]
            self.all_intensity_model_names = arrays[1]
            self.all_model_brightparams = arrays[2]
            self.all_model_names = arrays[3]
            self.total_models_count = arrays[4]
            fileloading.writeDocString(self.sub_paths["meta"] + "IntensityModelsGuide.txt",
                                       self.intensityModelDocString())
            fileloading.writeDocString(self.sub_paths["meta"] + "AllModelsGuide.txt",
                                       self.totalModelDocString())

        else:
            grid_types = {
                0: self.type0grid,
                1: self.type1Grid,
                2: self.type2Grid,
                3: self.type3Grid
            }

            self.all_intensity_model_names, self.all_intensity_model_brightparams = grid_types[self.run_type]()
            self.all_model_names = []
            self.all_model_brightparams = []
            self.total_models_count = 0

            fileloading.writeDocString(self.sub_paths["meta"] + "IntensityModelsGuide.txt",
                                       self.intensityModelDocString())

            # All model names
            self.creatAllModelNames()

            fileloading.writeDocString(self.sub_paths["meta"] + "AllModelsGuide.txt",
                                       self.totalModelDocString())

    def getSubPaths(self):
        return self.sub_paths

    def getModelBrightParams(self,model:str):

        return self.all_model_brightparams[list(self.all_model_names).index(model)]

    def type3Grid(self):
        all_brightparams = []  # list[tuple (name, bright parameters)]
        all_model_names = []
        for i in range(self.variable_param_ranges[self.var_intensity_grid_names[0]]):
            for j in range(self.variable_param_ranges[self.var_intensity_grid_names[1]]):
                for k in range(self.variable_param_ranges[self.var_intensity_grid_names[2]]):
                    current_model_brightparams = {}  # brightprams for current model

                    if i < 9:
                        place1 = str(i + 1)
                    else:
                        place1 = str(i + 1) + '-'

                    if j < 9:
                        place2 = str(j + 1)
                    elif j >= 9 and i >= 9:
                        place2 = str(j + 1) + '-'
                    else:
                        place2 = '-' + str(j + 1)
                    

                    if k < 9:
                        place3 = str(k + 1)
                    elif k >= 9 and j >= 9:
                        place3 = str(k + 1) + '-'
                    else:
                        place3 = '-' + str(k + 1)

                    current_model_name = "Model_" + place1 + place2 + place3

                    # Fill out the constant parameters
                    for key in list(self.constant_params):
                        current_model_brightparams[key] = self.constant_params[key]

                    current_model_brightparams[self.var_intensity_grid_names[0]] = \
                        self.intensity_grid_params[self.var_intensity_grid_names[0]][i]
                    current_model_brightparams[self.var_intensity_grid_names[1]] = \
                        self.intensity_grid_params[self.var_intensity_grid_names[1]][j]
                    current_model_brightparams[self.var_intensity_grid_names[2]] = \
                        self.intensity_grid_params[self.var_intensity_grid_names[2]][k]

                    all_brightparams += [current_model_brightparams]
                    all_model_names += [current_model_name]

        return all_model_names, all_brightparams

    def type2Grid(self):
        all_brightparams = []  # list[tuple (name, bright parameters)]
        all_model_names = []
        for i in range(self.variable_param_ranges[self.var_intensity_grid_names[0]]):
            for j in range(self.variable_param_ranges[self.var_intensity_grid_names[1]]):
                current_model_brightparams = {}  # brightprams for current model

                if i < 9:
                    place1 = str(i + 1)
                else:
                    place1 = str(i + 1) + '-'

                if j < 9:
                    place2 = str(j + 1)
                elif j >= 9 and i >= 9:
                    place2 = str(j + 1) + '-'
                else:
                    place2 = '-' + str(j + 1)
                current_model_name = "Model_" + place1 + place2

                # Fill out the constant parameters
                for key in list(self.constant_params):
                    current_model_brightparams[key] = self.constant_params[key]

                current_model_brightparams[self.var_intensity_grid_names[0]] = \
                    self.intensity_grid_params[self.var_intensity_grid_names[0]][i]
                current_model_brightparams[self.var_intensity_grid_names[1]] = \
                    self.intensity_grid_params[self.var_intensity_grid_names[1]][j]

                all_brightparams += [current_model_brightparams]
                all_model_names += [current_model_name]

        return all_model_names, all_brightparams

    def type1Grid(self):
        all_brightparams = []
        all_model_names = []
        for i in range(self.variable_param_ranges[self.var_intensity_grid_names[0]]):
            current_model_brightparams = {}  # dict[key] = for value of [key]/varied parameter

            current_model_name = "Model_" + str(i + 1)
            for key in list(self.constant_params):
                current_model_brightparams[key] = self.constant_params[key]

            current_model_brightparams[self.var_intensity_grid_names[0]] = \
                self.intensity_grid_params[self.var_intensity_grid_names[0]][i]

            all_brightparams += [current_model_brightparams]
            all_model_names += [current_model_name]

        return all_model_names, all_brightparams

    def type0grid(self):
        current_model_brightparams = {}  # dict[key] = for value of [key]/varied parameter
        current_model_name = "Model_" + str(1)
        for key in list(self.constant_params):
            current_model_brightparams[key] = self.constant_params[key]

        return [current_model_name], [current_model_brightparams]

    def createGeoGrid(self):
        for i in range(len(self.geo_grid_names)):
            print("\nCreating Geo Model '{}' \n".format(self.geo_grid_names[i]))
            print(line)
            # Create normalizing lensing band_____________________________________________________________________________________________

            normal_param_name = self.geo_grid_names[i] + "_Normalizing"

            fileloading.loadGeoModel(normal_param_name, self.run)
            importlib.reload(params)
            spin_case = params.spin_case
            i_case = params.i_case

            new_lband = self.sub_paths["GeoDoth5Path"] + normal_param_name + "Lensing" + ".h5"
            new_rtray = self.sub_paths["GeoDoth5Path"] + normal_param_name + "RayTracing" + ".h5"
            new_mag_angle = self.sub_paths["GeoDoth5Path"] + normal_param_name + "MagAng" + ".h5"
            print("Computation of {} Lensing Bands".format(normal_param_name))
            """Computation of the lensing bands______________________________________________________"""
            subprocess.run(['python3 ' + EZPaths.aartPath + '/lensingbands.py '], shell=True)
            fnbands = path + "LensingBands_a_%s_i_%s.h5" % (params.spin_case, params.i_case)

            '''Analytical Ray-tracing______________________________________________________'''
            print("Analytical Raytracing of {}".format(normal_param_name))
            subprocess.run(['python3 ' + EZPaths.aartPath + '/raytracing.py '], shell=True)
            fnrays1 = path + "Rays_a_%s_i_%s.h5" % (spin_case, i_case)
            
            """Computation of the Magnetic Angle______________________________________________________"""
            subprocess.run(['python3 ' + EZPaths.aartPath + '/magneticangle.py '], shell=True)
            fn= path + "MagneticAngle_a_%s_i_%s.h5"%(spin_case,i_case)

            # Move lensing bands and raytracing bands
            subprocess.run(["mv " + fnbands + ' ' + new_lband], shell=True)
            subprocess.run(["mv " + fnrays1 + ' ' + new_rtray], shell=True)
            subprocess.run(["mv " + fn + ' ' + new_mag_angle], shell=True)

            # Actual Model Params__________________________________________________________________________________________________________
            fileloading.loadGeoModel(self.geo_grid_names[i], self.run)
            importlib.reload(params)
            spin_case = params.spin_case
            i_case = params.i_case


            new_lband = self.sub_paths["GeoDoth5Path"] + self.geo_grid_names[i] + "Lensing" + ".h5"
            new_rtray = self.sub_paths["GeoDoth5Path"] + self.geo_grid_names[i] + "RayTracing" + ".h5"
            new_mag_angle = self.sub_paths["GeoDoth5Path"] + self.geo_grid_names[i] + "MagAng" + ".h5"
            print("Computation of {} Lensing Bands".format(self.geo_grid_names[i]))
            """Computation of the lensing bands______________________________________________________"""
            subprocess.run(['python3 ' + EZPaths.aartPath + '/lensingbands.py '], shell=True)
            fnbands = path + "LensingBands_a_%s_i_%s.h5" % (params.spin_case, params.i_case)

            '''Analytical Ray-tracing______________________________________________________'''
            print("Analytical Raytracing of {}".format(self.geo_grid_names[i]))
            subprocess.run(['python3 ' + EZPaths.aartPath + '/raytracing.py '], shell=True)
            fnrays1 = path + "Rays_a_%s_i_%s.h5" % (spin_case, i_case)


            """Computation of the Magnetic Angle______________________________________________________"""
            subprocess.run(['python3 ' + EZPaths.aartPath + '/magneticangle.py '], shell=True)
            fn= path + "MagneticAngle_a_%s_i_%s.h5"%(spin_case,i_case)


            # Move lensing bands and raytracing bands
            subprocess.run(["mv " + fnbands + ' ' + new_lband], shell=True)
            subprocess.run(["mv " + fnrays1 + ' ' + new_rtray], shell=True)
            subprocess.run(["mv " + fn + ' ' + new_mag_angle], shell=True)
            print("Files Moved and created: ")
            print("     ",new_lband)
            print("     ",new_rtray)
            print("     ",new_mag_angle)

    def creatAllModelNames(self):
        k = 0
        print("Creating all total model names")
        for j in range(len(self.geo_grid_names)):
            current_geo_model = self.geo_grid_names[j]
            fileloading.loadGeoModel(current_geo_model, self.run)

            for i in range(len(self.all_intensity_model_brightparams)):
                # String Names
                current_intent_name = self.all_intensity_model_names[i]
                self.all_model_names += [current_geo_model + current_intent_name.replace("Model", "")]

                newbp = self.all_intensity_model_brightparams[i].copy()
                self.all_model_brightparams += [newbp]
                self.total_models_count += 1
                k += 1

    def createBrightParamsDict(self):
        diction = {}
        for i in range(len(self.all_model_brightparams)):
            # String Names
            diction[self.all_model_names[i]] = self.all_model_brightparams[i]
            # ________________________________
        return diction

    def intensityModelDocString(self):
        string = line + line + line + "Total Number of Intensity Models: " + str(
            self.total_intensity_models_count) + '\n' + "Constant Params: " + '\n'
        for key in list(self.constant_params):
            string += breaker + key + ": " + str(self.constant_params[key]) + '\n'

        string += line
        for i in range(len(self.all_intensity_model_names)):
            string += line + self.all_intensity_model_names[i] + '\n'
            current_model = self.all_intensity_model_brightparams[i]
            for k in range(len(self.var_intensity_grid_names)):
                string += (breaker + self.var_intensity_grid_names[k] + ": "
                           + str(current_model[self.var_intensity_grid_names[k]]) + '\n')

        return string + line + line + line


    def totalModelDocString(self):

        intensity_model_string = line_small + line_small + line_small + "Total Number of Models: " + str(
            self.total_models_count) + '\n' + "Constant Params: " + '\n'
        for key in list(self.constant_params):
            if key != "n_th0":
                intensity_model_string += breaker + key + ": " + str(self.string_order[key].format(
                    self.constant_params[key])) + '\n'

        intensity_model_string += line_small
        for i in range(len(self.all_model_names)):
            current_name = self.all_model_names[i]
            current_model = self.all_model_brightparams[i]
            intensity_model_string += line_small + current_name + '\n'

            for k in range(len(self.var_intensity_grid_names)):
                intensity_model_string += (breaker + self.var_intensity_grid_names[k] + ": " +
                                           str(self.string_order[self.var_intensity_grid_names[k]].format(
                                               current_model[self.var_intensity_grid_names[k]]) + '\n'
                                               ))

            intensity_model_string += breaker + "n_th0: " + str(current_model["n_th0"]) + '\n'

        intensity_model_string += line_small + line_small + line_small

        geo_models_string = "\nGEOMODELS\n"
        for i in range(len(self.geo_grid_names)):
            geo_models_string += self.geo_grid_names[i] + "| "
            for j in range(len(self.geo_grid_params[i][0])):
                geo_models_string += '\n    ' + self.geo_grid_params[i][0][j] + ": " + self.geo_grid_params[i][1][j]
            geo_models_string += '\n'

        return intensity_model_string + geo_models_string

    """ Clean _______________________________________________________________________________________________________"""

    def creatIntensityGrid(self, action, do_list=None, isContinuous=False,frequency_list=None,average=True,analyze=False,keep_freq_list=None,blurr_list=None,blurr=False,kernal=20):

        all_230_total_jy_thin = []
        all_230_total_jy_thick = []

        k = 0
        print(line)
        print(line)
        print(line)
        print("Running Intensity Grid for " + self.run)
        print("All GeoGrid Names:  " + "\n" + str(self.geo_grid_names))
        for j in range(len(self.geo_grid_names)):
            current_geo_model = self.geo_grid_names[j]
            fileloading.loadGeoModel(current_geo_model, self.run)
            # lband = self.sub_paths["GeoDoth5Path"] + current_geo_model + "Lensing" + ".h5"
            # rtray = self.sub_paths["GeoDoth5Path"] + current_geo_model + "RayTracing" + ".h5"

            normlband = self.sub_paths["GeoDoth5Path"] + current_geo_model + "_Normalizing" + "Lensing" + ".h5"
            normrtray = self.sub_paths["GeoDoth5Path"] + current_geo_model + "_Normalizing" + "RayTracing" + ".h5"
            normMagAng = self.sub_paths["GeoDoth5Path"] + current_geo_model + "_Normalizing" + "MagAng" + ".h5"

            for i in range(len(self.all_intensity_model_brightparams)):
                # String Names
                current_total_name = self.all_model_names[k]
                # ________________________________
                current_bp = self.all_model_brightparams[k]

                parent_model_path = self.sub_paths["intensityPath"] + current_total_name + "/"
                current_model_file = parent_model_path + "clean/"

                fileloading.creatSubDirectory(parent_model_path,
                                              "for {} intensities".format(current_total_name), kill_policy=False)

                preform_model = fileloading.crossContinousDoAnalysis(
                    current_total_name, do_list, current_model_file, isContinuous)
                if preform_model:
                    print("Creating Intensity Movie for Model ", current_total_name)
                    print(long_line)

                    if do_list is None:
                        kill_policy = True
                    else:
                        kill_policy = False
                    fileloading.creatSubDirectory(current_model_file,
                                                  "for {} intensities".format(current_total_name),
                                                  kill_policy=kill_policy)

                    if not self.already_normalized_brightparams:
                        print("\n" + "Normalizing " + current_total_name + "\n")
                        print(long_line)

                        current_bp["n_th0"] = normalizingBrightparams.normalize(normlband, normrtray,normMagAng,current_bp,self.funckeys)

                        print("\n" + current_total_name + " normalized with a value of n_th0="
                              + str(current_bp["n_th0"]) + "\n")
                        self.all_model_brightparams[k] = current_bp

                    if self.run_type == 0:
                        run_type_arg = 1
                    else:
                        run_type_arg = self.run_type
                    intermodel_data = movieMakerIntensity.intensity_movie(
                        action, self.sub_paths, current_total_name, run_type_arg, current_bp,frequency_list,self.funckeys)

                    print(
                        "\nTotal flux at 230GHz for Optically Thin Assumption: " + str(
                            intermodel_data["thin_total_flux"]))
                    print("Total flux at 230GHz for Full Solution: " + str(intermodel_data["thick_total_flux"]) + "\n")
                    if do_list is None:
                        all_230_total_jy_thin += [intermodel_data["thin_total_flux"]]
                        all_230_total_jy_thick += [intermodel_data["thick_total_flux"]]
                    else:
                        print("Not all Models are being analyzed, skipping collection of intermodel data")

                    # Analysis_____________
                    if analyze:
                        if keep_freq_list is None:
                            keep_freq_list = np.arange(action["start"],action["stop"] - action["step"],action["step"])
                        #     print(keep_freq_list)
                        # str_keep_freq_list = []
                        # for i in range(len(keep_freq_list)):
                        #     str_keep_freq_list += ["{:.5e}".format(keep_freq_list[i])]

                        self.analyzeModel(current_total_name,action,current_bp,average)

                        if blurr:
                            if blurr_list is not None:
                                if current_total_name in blurr_list:
                                    do_blurr = True
                                else:
                                    do_blurr = False
                            else:
                                do_blurr = True

                            if do_blurr:
                                if type(kernal) is not list:
                                    kernal = [kernal]
                                for i in range(len(kernal)):
                                    current_kernal = kernal[i]
                                    movieMakerIntensity.blurr_intensity_movie(
                                        action, self.sub_paths, current_total_name, run_type_arg, current_bp,None,current_kernal,self.funckeys
                                    )
                                    
                                    movieMakerIntensity.blurrImageAnalysis(
                                        action, self.sub_paths, current_total_name, current_bp, None, current_kernal
                                    )
                                    blurr_dir = self.sub_paths["intensityPath"] + current_total_name + "/" + "blurr" + str(current_kernal) + "/"
                                    self.deleteNu0Points(blurr_dir,keep_freq_list)
                        
                        self.deleteNu0Points(self.sub_paths["intensityPath"] + current_total_name + "/" + "clean/", keep_freq_list)
                else:
                    print("{} marked for skipping...".format(current_total_name))
                k += 1
        if do_list is None:
            if not self.already_normalized_brightparams:
                print("Saving Normalized Parameters")
                file_paths = [
                    self.sub_paths["runWideNumpy"] + "all_intensity_model_brightparams",
                    self.sub_paths["runWideNumpy"] + "all_intensity_model_names",
                    self.sub_paths["runWideNumpy"] + "all_model_brightparams",
                    self.sub_paths["runWideNumpy"] + "all_model_names",
                    self.sub_paths["runWideNumpy"] + "total_models_count"
                ]
                arrays = [
                    self.all_intensity_model_brightparams,
                    self.all_intensity_model_names,
                    self.all_model_brightparams,
                    self.all_model_names,
                    self.total_models_count
                ]
                k = 0
                for file in file_paths:
                    np.save(file, arrays[k])
                    k += 1

                fileloading.writeDocString(self.sub_paths["meta"] + "IntensityModelsGuide.txt",
                                        self.intensityModelDocString())
                fileloading.writeDocString(self.sub_paths["meta"] + "AllModelsGuide.txt",
                                        self.totalModelDocString())

            self.already_normalized_brightparams = True
            # Numpy saving________________________________________________

            all_230_total_jy_thin_numpy_name = self.sub_paths["meta"] + "thin_total_flux"
            all_230_total_jy_thick_numpy_name = self.sub_paths["meta"] + "thick_total_flux"

            np.save(all_230_total_jy_thin_numpy_name, np.array(all_230_total_jy_thin))
            np.save(all_230_total_jy_thick_numpy_name, np.array(all_230_total_jy_thick))
        else:
            print("Not all Models are being analyzed, skipping saving of intermodel data")

    def creatRadialProfiles(self, action, do_list=None, isContinuous=False, frequency_list=None):

        k = 0
        print(line)
        print(line)
        print(line)
        print("Creating radial profiles for " + self.run)

        for i in range(len(self.all_model_brightparams)):
            # String Names
            current_total_name = self.all_model_names[i]
            # ________________________________
            current_bp = self.all_model_brightparams[i]

            # File Creation_____________________________
            parent_model_path = self.sub_paths["intensityPath"] + current_total_name + "/"
            current_model_file = parent_model_path + "clean/"
            final_data_path = current_model_file + "numpy/"

            preform_model = fileloading.crossContinousDoAnalysis(
                current_total_name, do_list, final_data_path, isContinuous)

            if preform_model:
                fileloading.creatSubDirectory(final_data_path,
                                              "final image path for {}".format(current_total_name), kill_policy=False)

                # __________________________________________

                print("Analyzing Intensity Movie for Model ", current_total_name)
                print(long_line)
                movieMakerIntensity.profileAnalysis(action,
                                                    self.sub_paths,
                                                    current_total_name,
                                                    current_bp,
                                                    frequency_list)
            else:
                print(current_total_name + " marked for skipping...")

    def intensityGridAnalysis(self, action, do_list=None, isContinuous=False, average=True):
        print(line)
        print(line)
        print(line)
        print("Analyzing Intensity Grid for " + self.run)
        for i in range(len(self.all_model_brightparams)):
            # String Names
            current_total_name = self.all_model_names[i]
            # ________________________________
            current_bp = self.all_model_brightparams[i]

            parent_model_path = self.sub_paths["intensityPath"] + current_total_name + "/"
            current_model_file = parent_model_path + "clean/"
            preform_model = fileloading.crossContinousDoAnalysis(
                current_total_name, do_list, current_model_file, isContinuous)

            if preform_model:
                self.analyzeModel(current_total_name,action,current_bp,average)

                # print("Analyzing Intensity Movie for Model ", current_total_name)
                # print(long_line)
                # movieMakerIntensity.imageAnalysis(
                #     action, self.sub_paths, current_total_name, current_bp,average=average
                # )
            else:
                print(current_total_name + " marked for skipping...")
    
    def analyzeModel(self,model_name,action,bp,average):
        print("Analyzing Intensity Movie for Model ", model_name)
        print(long_line)
        movieMakerIntensity.imageAnalysis(
            action, self.sub_paths, model_name, bp,average=average
        )

    def deleteNu0Points(self,dire,keep_list):
        files = fileloading.filesInDir(dire)

        str_keep_list = []
        for i in range(len(keep_list)):
            str_keep_list += ["{:.5e}".format(keep_list[i])]

            
        for file in files:
            full_file = dire + file
            if not fileloading.strContains(file,str_keep_list):
                os.remove(full_file)

    def interModelAnalysis(self, action, do_list=None, average=True):
        k = 0
        print(line)
        print(line)
        print(line)
        print("Running Intermodel Analysis for " + self.run)
        print("All GeoGrid Names:  " + "\n" + str(self.geo_grid_names))

        # Doc String_______________________________________
        string = line_small + line_small + line_small + "Total Number of Models: " + str(
            self.total_models_count) + '\n' + "Constant Params: " + '\n'
        for key in list(self.constant_params):
            if key != "n_th0":
                if (key == "mass") or (key == "nu0"):
                    string += breaker + key + ": " + str("{:.5e}".format(self.constant_params[key])) + '\n'
                elif key == "theta_b":
                    string += breaker + key + ": " + str("{:.5f}".format(self.constant_params[key])) + '\n'
                else:
                    string += breaker + key + ": " + str(self.constant_params[key]) + '\n'

        string += line_small

        # _______________________________________

        for model in self.all_model_names:

            if do_list is not None:
                preform_model = fileloading.doListAnalysis(model,do_list)
            else:
                preform_model = True

            if preform_model:
                print(line)
                print("Printing Analysis of " + model)
                # File Creation

                # if do_list is None:
                #     kill_policy = True
                # else:
                #     kill_policy = False

                if self.run_type == 0:
                    amount_to_subtract = 1
                else:
                    amount_to_subtract = self.run_type

                current_geo_model = model.split("_")[0]
                fileloading.loadGeoModel(current_geo_model, self.run)
                # lband = self.sub_paths["GeoDoth5Path"] + current_geo_model + "Lensing" + ".h5"
                # rtray = self.sub_paths["GeoDoth5Path"] + current_geo_model + "RayTracing" + ".h5"

                # Construct Shadows___________________________________________________________________
                # a = params.spin_case
                # inc = params.i_case * np.pi / 180  # inclination angle
                # rh = 1 + np.sqrt(1 - a ** 2)  # event horizon
                # # angles to sample
                # varphis = np.linspace(-180, 179, 360) * np.pi / 180
                #
                # # generate inner shadow (n=0) curve with kgeo
                # data_inner = kgeo.equatorial_lensing.rho_of_req(a, inc, rh, mbar=0, varphis=varphis)
                # (_, rhos_inner, alphas_inner, betas_inner) = data_inner
                #
                # r_inner = image_tools.curve_params(varphis, rhos_inner)
                #
                # # generate outer shadow (n=inf) curve with kgeo
                # data_outer = kgeo.equatorial_lensing.rho_of_req(a, inc, rh, mbar=5, varphis=varphis)
                # (_, rhos_outer, alphas_outer, betas_outer) = data_outer
                #
                # r_outer = image_tools.curve_params(varphis, rhos_outer)
                # ___________________________________________________________________

                '''Data Reading----------------------------------'''
                data_path = self.sub_paths["intensityPath"] + model + "/clean/numpy/"

                print("Reading: " + data_path)

                if not average:
                    data_path += "FalseAvg_"

                x_variable = np.load(data_path + "x_variable.npy")
                janksys_thick = np.load(data_path + "janksys_thick.npy")
                janksys_thin = np.load(data_path + "janksys_thin.npy")
                # mean_radii_Thin = np.load(data_path + "mean_radii_Thin.npy")
                mean_radii_Thick = np.load(data_path + "mean_radii_Thick.npy")
                # radii_I0_Thin = np.load(data_path + "radii_I0_Thin.npy")
                # radii_I1_Thin = np.load(data_path + "radii_I1_Thin.npy")
                # radii_I2_Thin = np.load(data_path + "radii_I2_Thin.npy")
                # radii_Full_Thin = np.load(data_path + "radii_Full_Thin.npy")
                # radii_FullAbsorption_Thick = np.load(data_path + "radii_FullAbsorption_Thick.npy")
                # radii_I0_Thick = np.load(data_path + "radii_I0_Thick.npy")
                # radii_I1_Thick = np.load(data_path + "radii_I1_Thick.npy")
                # radii_I2_Thick = np.load(data_path + "radii_I2_Thick.npy")
                # theta = np.load(data_path + "theta.npy")
                # mean_optical_depth_I0 = np.load(data_path + "mean_optical_depth_I0.npy")
                # mean_optical_depth_I1 = np.load(data_path + "mean_optical_depth_I1.npy")
                # mean_optical_depth_I2 = np.load(data_path + "mean_optical_depth_I2.npy")

                num_of_intensity_points = janksys_thin[:, 0].shape[0]
                print("Number of Intensity Points: ", num_of_intensity_points)

                xaxis = np.array(x_variable) / scale_label[action['var']]
                # one_M = ilp.rg_func(brightparams["mass"] * u.g).to(u.m)
                # M2uas = np.arctan(one_M.value / dBH) / muas_to_rad

                # Points of Interest

                conv_1 = ilp.ring_convergance(xaxis,mean_radii_Thick[:, 2], mean_radii_Thick[:, 3], 3)

                flux_peak_thin = ilp.function_peak(xaxis,janksys_thin[:, 3])

                flux_peak_thick = ilp.function_peak(xaxis,janksys_thick[:, 3])

                # String Data
                string += line_small + model + '\n'

                string += (breaker + R"Full RTE Solution $\nu_{peak}$" + ": " + str("{:.5f}".format(flux_peak_thick)) + 'GHz\n')
                # string += (breaker + R"Optically Thin Assumption $\nu_{peak}$" + ": " + str(flux_peak_thin) + '\n')
                string += (breaker + R"Full RTE Solution $\nu_{conv}$" + ": " + str("{:.5f}".format(conv_1)) + 'GHz\n')

        fileloading.writeDocString(self.sub_paths["meta"] + "interModelAnalysis.txt",
                                   string)

    def graphCreation(self, action, do_list=None, isContinuous=False,average=True,doFullImages=True):
        """

        Args:
            action: = {"var":str, "start":float, "stop":float, "step":float}

        Returns:

        """
        # PLOT SETTINGS__________________________________________________
        plt.rcParams.update({
            'font.size': 14,  # Set font size to 11pt
            'axes.labelsize': 14,  # -> axis labels
            'legend.fontsize': 12,  # -> legends
            'text.usetex': True,
            'text.latex.preamble': (  # LaTeX preamble
                r'\usepackage{lmodern}'
            ),
            'font.family': 'Latin Modern Roman',
        })

        # ________________________________________________________________
        print(line)
        print(line)
        print("Initializing Graph Creation")
        thin_total_flux = np.load(self.sub_paths["meta"] + "thin_total_flux.npy")
        thick_total_flux = np.load(self.sub_paths["meta"] + "thick_total_flux.npy")
        # current_intensity_model_name = intensity_models[i][0]
        # current_intensity_model = intensity_models[i][1]
        # full_current_model_name = current_geo_model + current_intensity_model_name.replace("Model", "")
        j = 0
        if do_list is None:
            hist_flux_peaks_thins = []
            hist_flux_peaks_thicks = []
            hist_convs = []
        dim = [10, 8]
        for model in self.all_model_names:

            fluxVNu_path = self.sub_paths["fluxPath"] + model
            radVNu_path = self.sub_paths["radPath"] + model
            image_path = self.sub_paths["imagePath"] + model
            Optical_depth_path = self.sub_paths["opticalDepth"] + model
            radVVarphi_path = self.sub_paths["RadVVarphi"] + model
            fluxVRadii_path = self.sub_paths["fluxVRadii"] + model

            file_creation = [
                fluxVNu_path,
                radVNu_path,
                image_path,
                Optical_depth_path,
                radVVarphi_path,
                fluxVRadii_path
            ]

            preform_model = fileloading.crossContinousDoAnalysis(
                model, do_list, fluxVNu_path, isContinuous)

            if preform_model:
                print(line)
                print("Graphing " + model)
                # File Creation
                for i in range(len(file_creation)):
                    fileloading.creatSubDirectory(file_creation[i], kill_policy=False)

                fluxVNu_path +="/clean/"
                radVNu_path += "/clean/"
                image_path += "/clean/"
                Optical_depth_path += "/clean/"
                radVVarphi_path += "/clean/"
                fluxVRadii_path += "/clean/"

                file_creation = [
                    fluxVNu_path,
                    radVNu_path,
                    image_path,
                    Optical_depth_path,
                    radVVarphi_path,
                    fluxVRadii_path
                ]

                if do_list is None:
                    kill_policy = True
                else:
                    kill_policy = False

                for i in range(len(file_creation)):
                    if (file_creation[i] == image_path) and doFullImages is False:
                        fileloading.creatSubDirectory(file_creation[i], kill_policy=False)
                    else:
                        fileloading.creatSubDirectory(file_creation[i], kill_policy=kill_policy)

                fluxVNu_path += "Clean_"
                radVNu_path += "Clean_"
                image_path += "Clean_"
                Optical_depth_path += "Clean_"
                radVVarphi_path += "Clean_"
                fluxVRadii_path += "Clean_"

                if not average:
                    fluxVNu_path += "FalseAvg_"
                    radVNu_path += "FalseAvg_"
                    image_path += "FalseAvg_"
                    Optical_depth_path += "FalseAvg_"
                    radVVarphi_path += "FalseAvg_"
                    fluxVRadii_path += "FalseAvg_"

                if self.run_type == 0:
                    amount_to_subtract = 1
                else:
                    amount_to_subtract = self.run_type

                current_geo_model = model.split("_")[0]
                fileloading.loadGeoModel(current_geo_model, self.run)
                lband = self.sub_paths["GeoDoth5Path"] + current_geo_model + "Lensing" + ".h5"
                rtray = self.sub_paths["GeoDoth5Path"] + current_geo_model + "RayTracing" + ".h5"

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
                data_path = self.sub_paths["intensityPath"] + model + "/clean/numpy/"

                if not average:
                    data_path += "FalseAvg_"

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

                xaxis = np.array(x_variable) / scale_label[action['var']]
                # one_M = ilp.rg_func(brightparams["mass"] * u.g).to(u.m)
                # M2uas = np.arctan(one_M.value / dBH) / muas_to_rad

                # Points of Interest

                conv_1 = ilp.ring_convergance(xaxis,mean_radii_Thick[:, 2], mean_radii_Thick[:, 3], 3)
                flux_peak_thin = ilp.function_peak(xaxis,janksys_thin[:, 3])
                flux_peak_thick = ilp.function_peak(xaxis,janksys_thick[:, 3])

                if do_list is None:
                    hist_flux_peaks_thins += [flux_peak_thin]
                    hist_flux_peaks_thicks += [flux_peak_thick]
                    hist_convs += [conv_1]

                poi = {
                    "r_outer": r_outer,
                    "flux_peak_thin": flux_peak_thin,
                    "flux_peak_thick": flux_peak_thick,
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
                                          [mean_optical_depth_I0, mean_optical_depth_I1, mean_optical_depth_I2],
                                          poi, conv_1_style, flux_peak_style, action)

                figname = Optical_depth_path + model + "OpticalDepth.jpeg"
                plt.savefig(figname, bbox_inches='tight')
                print("Image '{}' Created".format(figname))
                plt.close()

                '''Full Images----------------------------------'''
                if doFullImages:
                    parent_model_path = self.sub_paths["intensityPath"] + model + "/"
                    current_model_file = parent_model_path + "clean/"
                    k = action["start"]
                    print("Constructing Full images for " + model)
                    for i in range(num_of_intensity_points):
                        brightparams = self.all_model_brightparams[j]
                        brightparams["nu0"] = k
                        print("Full image production for intensity frame: ", i)
                        print(R"Observation frequency $\nu=$", k)

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
                        rmax = I0.shape[0] * .4
                        thin_intensity = [I0, I1, I2, I0 + I1 + I2]
                        thick_intensity = [I0_Absorb, I1_Absorb, I2_Absorb, Absorbtion_Image]
                        thin_radii = [radii_I0_Thin[i, :], radii_I1_Thin[i, :], radii_I2_Thin[i, :],
                                      radii_Full_Thin[i, :]]
                        thick_radii = [radii_I0_Thick[i, :], radii_I1_Thick[i, :],
                                       radii_I2_Thick[i, :], radii_FullAbsorption_Thick[i, :]]

                        fig, (ax0, ax1) = plt.subplots(1, 2, figsize=[15, 7], dpi=400)

                        astroPloting.fullImage(fig, ax0, ax1, lim0, thin_intensity, thick_intensity,
                                               thin_radii, thick_radii, theta)

                        ax0.text(-9, 8.5, var_label[action["var"]]
                                 + str(round(x_variable[i] / scale_label[action["var"]], 2))
                                 + ' ' + units_label[action["var"]], fontsize=12, color="w")

                        pltname = (image_path + 'FullImage_' + str(i) + "_Nu_"
                                   + str(round(x_variable[i] / scale_label[action["var"]], 2)) + ".jpeg")
                        plt.savefig(pltname, bbox_inches='tight')
                        print("Jpeg Created:  " + pltname)
                        plt.close()

                        # Get total jansky
                        # VARPHI__________________________________________________-

                        '''IntensityVSRadii________________________________________________________________'''
                        fig, dum = plt.subplots(1, 2, figsize=[15, 7], dpi=400)
                        ax0 = plt.subplot(1, 2, 1)
                        ax1 = plt.subplot(1, 2, 2)

                        astroPloting.IntensityVSRadii(fig, ax0, ax1, params.limits, thin_intensity,
                                                      thick_intensity, rmax)
                        #
                        # ax1.text(2, 1.01, var_label[action["var"]]
                        #          + str(round(x_variable[i] / scale_label[action["var"]], 2))
                        #          + ' ' + units_label[action["var"]], fontsize=12, color="k")

                        pltname = (fluxVRadii_path + 'IntVRad_' + str(i) + "_Nu_"
                                   + str(round(k / scale_label[action["var"]], 2)) + ".jpeg")
                        plt.savefig(pltname, bbox_inches='tight')
                        print("Image '{}' Created".format(pltname))
                        plt.close()

                        '''RadVSVarphiType2________________________________________________________________'''
                        fig, dum = plt.subplots(1, 2, figsize=[15, 7], dpi=400)
                        ax0 = plt.subplot(1, 2, 1)
                        ax1 = plt.subplot(1, 2, 2)

                        astroPloting.radiiVSVarphi(fig, ax0, ax1, params.limits, thin_intensity,average=average)

                        # ax0.text(2, 1.01, var_label[action["var"]]
                        #          + str(round(x_variable[i] / scale_label[action["var"]], 2))
                        #          + ' ' + units_label[action["var"]], fontsize=12, color="k")

                        pltname = (radVVarphi_path + 'radVVarphu_' + str(i) + "_Nu_"
                                   + str(round(k / scale_label[action["var"]], 2)) + ".jpeg")
                        plt.savefig(pltname, bbox_inches='tight')
                        print("Image '{}' Created".format(pltname))
                        plt.close()

                        k += action['step']
            else:
                print(model + " marked for skipping...")

            j += 1  # marker for which brightparams to use

        # histograms
        if do_list is None:
            print(line)
            print("Creating Histograms")
            peak_hist_thin_path = self.sub_paths["peakHistThin"] + "clean"
            peak_hist_thick_path = self.sub_paths["peakHistThick"] + "clean"
            conv_hist_path = self.sub_paths["convHist"] + "clean"
            total_flux_path = self.sub_paths["totalFlux"] + "clean"

            bar_xaxis = np.arange(len(self.all_model_names))
            bar_labels = self.all_model_names
            for i in range(len(bar_xaxis)):
                bar_labels[i] = bar_labels[i].replace("Model", "")
            """Flux Peaks_____________________________________________________________________"""
            # Thin_________________

            fig, ax = plt.subplots(1, 1, figsize=dim, dpi=400)

            astroPloting.histogram(ax, hist_flux_peaks_thins, "Flux Peak location (GHz)",
                                   "Optically Thin Assumption Frequency")

            figname = peak_hist_thin_path + "FluxPeakThin.jpeg"
            plt.savefig(figname, bbox_inches='tight')
            print("Image '{}' Created".format(figname))
            plt.close()

            # Thick_______________

            fig, ax = plt.subplots(1, 1, figsize=dim, dpi=400)

            astroPloting.histogram(ax, hist_flux_peaks_thicks, "Flux Peak location (GHz)", "Full Solution Frequency")

            figname = peak_hist_thick_path + "FluxPeakThick.jpeg"
            plt.savefig(figname, bbox_inches='tight')
            print("Image '{}' Created".format(figname))
            plt.close()

            # Thin Bar__________________
            fig, ax = plt.subplots(1, 1, figsize=dim, dpi=400)

            astroPloting.bar(ax, bar_xaxis, hist_flux_peaks_thins, "Optical Thin Assumption Peak Flux per model",
                             "Observation Frequency (GHz)", bar_labels)

            figname = peak_hist_thin_path + "FluxPeakPerThinModel.jpeg"
            plt.savefig(figname, bbox_inches='tight')
            print("Image '{}' Created".format(figname))
            plt.close()

            # Thick Bar_______________

            fig, ax = plt.subplots(1, 1, figsize=dim, dpi=400)

            astroPloting.bar(ax, bar_xaxis, hist_flux_peaks_thicks, "Optical Thin Assumption Peak Flux per model",
                             "Observation Frequency (GHz)", bar_labels)

            figname = peak_hist_thick_path + "FluxPeakPerThickModel.jpeg"
            plt.savefig(figname, bbox_inches='tight')
            print("Image '{}' Created".format(figname))
            plt.close()
            # Thick Bar

            """Conv_____________________________________________________________________"""
            # Conv Hist ______________________
            fig, ax = plt.subplots(1, 1, figsize=dim, dpi=400)

            astroPloting.histogram(ax, hist_convs, "Flux Peak location (GHz)", "Full Solution Frequency")

            figname = conv_hist_path + "convHistFreqPerObservationFreq.jpeg"
            plt.savefig(figname, bbox_inches='tight')
            print("Image '{}' Created".format(figname))
            plt.close()

            # Conv Bar______________________
            fig, ax = plt.subplots(1, 1, figsize=dim, dpi=400)

            astroPloting.bar(ax, bar_xaxis, hist_convs, "Convergence for cumulative on I2",
                             "Observation Frequency (GHz)", bar_labels)
            figname = conv_hist_path + "convHistPerModel.jpeg"
            plt.savefig(figname, bbox_inches='tight')
            print("Image '{}' Created".format(figname))
            plt.close()

            """230 total flux Thin___________________________________________________________"""
            fig, ax = plt.subplots(1, 1, figsize=dim, dpi=400)

            astroPloting.bar(ax, bar_xaxis, thin_total_flux, "Total Flux at 230GHz",
                             "Optically Thin Assumption Frequency", bar_labels)

            figname = total_flux_path + "thin.jpeg"
            plt.savefig(figname, bbox_inches='tight')
            print("Image '{}' Created".format(figname))
            plt.close()

            """230 total flux Thick___________________________________________________________"""
            fig, ax = plt.subplots(1, 1, figsize=dim, dpi=400)

            astroPloting.bar(ax, bar_xaxis, thick_total_flux, "Total Flux at 230GHz",
                             "Full Solution Frequency", bar_labels)

            figname = total_flux_path + "thick.jpeg"
            plt.savefig(figname, bbox_inches='tight')
            print("Image '{}' Created".format(figname))
            plt.close()
        else:
            print("Not all models were computed, skipping Histogram creation...")


    """ Blurr _______________________________________________________________________________________________________"""
    def blurrIntensityGrid(self, action, do_list=None, isContinuous=False,
                           blurr_frequency_list: list = None, blur_kernal: int = 20):

        all_230_total_jy_thin = []
        all_230_total_jy_thick = []

        print(line)
        print(line)
        print(line)
        print("Bluring Intensity Grid for " + self.run)
        for i in range(len(self.all_model_brightparams)):

            # String Names
            # current_intent_name = self.all_intensity_model_names[i]
            # self.all_model_names += [current_geo_model + current_intent_name.replace("Model", "")]
            current_total_name = self.all_model_names[i]
            # ________________________________
            current_bp = self.all_model_brightparams[i]

            if self.run_type == 0:
                run_type_arg = 1
            else:
                run_type_arg = self.run_type

            parent_model_path = self.sub_paths["intensityPath"] + current_total_name + "/"
            current_model_file = parent_model_path + "blurr" + str(blur_kernal) + "/"

            preform_model = fileloading.crossContinousDoAnalysis(
                current_total_name, do_list, current_model_file, isContinuous)

            if preform_model:
                print("Blurring " + current_total_name)
                intermodel_data = movieMakerIntensity.blurr_intensity_movie(
                    action, self.sub_paths, current_total_name, run_type_arg, current_bp,
                    blurr_frequency_list, blur_kernal,self.funckeys)

                print(
                    "\nTotal flux at 230GHz for Optically Thin Assumption: " + str(intermodel_data["thin_total_flux"]))
                print("Total flux at 230GHz for Full Solution: " + str(intermodel_data["thick_total_flux"]) + "\n")
                if do_list is None:
                    all_230_total_jy_thin += [intermodel_data["thin_total_flux"]]
                    all_230_total_jy_thick += [intermodel_data["thick_total_flux"]]
                else:
                    print("Not all Models are being analyzed, skipping collection of intermodel data")
            else:
                print(current_total_name + " marked for skipping...")
        # Numpy saving________________________________________________
        if do_list is None:
            all_230_total_jy_thin_numpy_name = self.sub_paths["meta"] + "blurr_thin_total_flux"
            all_230_total_jy_thick_numpy_name = self.sub_paths["meta"] + "blurr_thick_total_flux"

            np.save(all_230_total_jy_thin_numpy_name, np.array(all_230_total_jy_thin))
            np.save(all_230_total_jy_thick_numpy_name, np.array(all_230_total_jy_thick))
        else:
            print("Not all Models are being analyzed, skipping saving of intermodel data")

    def blurrIntensityGridAnalysis(self, action, do_list=None, isContinuous=False, blurr_frequency_list=None,
                                   blur_kernal=20,keep_freq_list=None):
        print(line)
        print(line)
        print(line)
        print("Analyzing blurred intensity grid for " + self.run)

        if keep_freq_list is None:    
            keep_freq_list = np.arange(action["start"],action["stop"] - action["step"],action["step"])
        for i in range(len(self.all_model_brightparams)):

            # String Names
            current_total_name = self.all_model_names[i]
            # ________________________________
            current_bp = self.all_model_brightparams[i]

            parent_model_path = self.sub_paths["intensityPath"] + current_total_name + "/"
            current_model_file =  parent_model_path + "blurr" + str(blur_kernal) + "/"

            preform_model = fileloading.crossContinousDoAnalysis(
                current_total_name, do_list, current_model_file, isContinuous)

            if preform_model:
                print("Analyzing blurred intensity movie for model ", current_total_name)

                movieMakerIntensity.blurrImageAnalysis(
                    action, self.sub_paths, current_total_name, current_bp, blurr_frequency_list, blur_kernal
                )
            else:
                print(current_total_name + " marked for skipping...")
            


    def blurrGraphCreation(self, action, do_list=None, isContinuous=False, blurr_frequency_list=None,
                           blur_kernal=20):
        """

        Args:
            action:  = {"var":str, "start":float, "stop":float, "step":float}
            do_list:
            isContinuous:
            blurr_frequency_list:
            blur_kernal:

        Returns:

        """
        # PLOT SETTINGS__________________________________________________
        plt.rcParams.update({
            'font.size': 14,  # Set font size to 11pt
            'axes.labelsize': 14,  # -> axis labels
            'legend.fontsize': 12,  # -> legends
            'text.usetex': True,
            'text.latex.preamble': (  # LaTeX preamble
                r'\usepackage{lmodern}'
            ),
            'font.family': 'Latin Modern Roman',
        })
        # ________________________________________________________________

        print(line)
        print(line)
        print("Initializing Blurr Graph Creation")
        thin_total_flux = np.load(self.sub_paths["meta"] + "thin_total_flux.npy")
        thick_total_flux = np.load(self.sub_paths["meta"] + "thick_total_flux.npy")

        j = 0
        hist_flux_peaks_thins = []
        hist_flux_peaks_thicks = []
        hist_convs = []
        dim = [10, 8]
        for model in self.all_model_names:

            print(line)
            print("Running " + model)

            # File Analysis__________________________________________________________________
            if self.run_type == 0:
                amount_to_subtract = 1
            else:
                amount_to_subtract = self.run_type

            current_geo_model = model.split("_")[0]
            fileloading.loadGeoModel(current_geo_model, self.run)

            fluxVNu_path = self.sub_paths["fluxPath"] + model
            radVNu_path = self.sub_paths["radPath"] + model
            image_path = self.sub_paths["imagePath"] + model
            Optical_depth_path = self.sub_paths["opticalDepth"] + model
            fluxVRadii_path = self.sub_paths["fluxVRadii"] + model
            radVVarphi_path = self.sub_paths["RadVVarphi"] + model

            file_creation = [fluxVNu_path, radVNu_path, image_path, Optical_depth_path]

            preform_model = fileloading.crossContinousDoAnalysis(
                model, do_list, fluxVNu_path, isContinuous)

            if preform_model:
                print(" Creating Graphs for " + model)

                for i in range(len(file_creation)):
                    fileloading.creatSubDirectory(file_creation[i], kill_policy=False)

                fluxVNu_path        += "/blurr" + str(blur_kernal) + "/"
                radVNu_path         += "/blurr" + str(blur_kernal) + "/"
                image_path          += "/blurr" + str(blur_kernal) + "/"
                Optical_depth_path  += "/blurr" + str(blur_kernal) + "/"
                fluxVRadii_path     += "/blurr" + str(blur_kernal) + "/"
                radVVarphi_path     += "/blurr" + str(blur_kernal) + "/"

                file_creation = [fluxVNu_path, radVNu_path, image_path, Optical_depth_path,fluxVRadii_path,radVVarphi_path]

                if do_list is None:
                    kill_policy = True
                else:
                    kill_policy = False

                for i in range(len(file_creation)):
                    fileloading.creatSubDirectory(file_creation[i], kill_policy=kill_policy)

                fluxVNu_path        += "Blurr" + str(blur_kernal) + "_"
                radVNu_path         += "Blurr" + str(blur_kernal) + "_"
                image_path          += "Blurr" + str(blur_kernal) + "_"
                Optical_depth_path  += "Blurr" + str(blur_kernal) + "_"
                fluxVRadii_path     += "Blurr" + str(blur_kernal) + "_"
                radVVarphi_path     += "Blurr" + str(blur_kernal) + "_"

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
                data_path = self.sub_paths["intensityPath"] + model + "/blurr" + str(blur_kernal) + "/"+ "numpy/"

                x_variable = np.load(data_path + "blurr_x_variable.npy")
                janksys_thick = np.load(data_path + "blurr_janksys_thick.npy")
                janksys_thin = np.load(data_path + "blurr_janksys_thin.npy")
                mean_radii_Thin = np.load(data_path + "blurr_mean_radii_Thin.npy")
                mean_radii_Thick = np.load(data_path + "blurr_mean_radii_Thick.npy")
                radii_I0_Thin = np.load(data_path + "blurr_radii_I0_Thin.npy")
                radii_FullAbsorption_Thick = np.load(data_path + "blurr_radii_FullAbsorption_Thick.npy")
                theta = np.load(data_path + "blurr_theta.npy")

                num_of_action_points = int((action["stop"] - action["start"]) / action["step"])
                print("Number of Intensity Points: ", num_of_action_points)

                xaxis = np.array(x_variable) / scale_label[action['var']]
                # one_M = ilp.rg_func(brightparams["mass"] * u.g).to(u.m)
                # M2uas = np.arctan(one_M.value / dBH) / muas_to_rad

                # Points of Interest

                conv_1 = ilp.ring_convergance(xaxis,mean_radii_Thick[:, 0],mean_radii_Thick[:, 0],3)

                flux_peak_thin = ilp.function_peak(xaxis, janksys_thin[:, 0])
                flux_peak_thick = ilp.function_peak(xaxis, janksys_thick[:, 0])

                hist_flux_peaks_thins += [flux_peak_thin]
                hist_flux_peaks_thicks += [flux_peak_thick]
                hist_convs += [conv_1]

                poi = {
                    "r_outer": r_outer,
                    "flux_peak_thin": flux_peak_thin,
                    "flux_peak_thick": flux_peak_thick,
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

                # _______________________________________________
                # _______________________________________________
                # ________________________________
                '''JANKSKY PLOTS----------------------------------'''

                fig, (ax, ax1) = plt.subplots(2, 1, figsize=dim, dpi=400, sharex=True,sharey=True)

                astroPloting.fluxThickThin(ax, ax1, xaxis, janksys_thin, janksys_thick,
                                           poi, conv_1_style, r_outer_style, flux_peak_style, action, blurr_policy=True)

                figname = fluxVNu_path + model + "_Flux.jpg"
                plt.savefig(figname, bbox_inches='tight')
                plt.close()
                print("Image '{}' Created".format(figname))

                # ______________________________________________
                # ______________________________________________
                # ___________________________________
                '''RADII PLOTS----------------------------------'''
                fig, (ax, ax1) = plt.subplots(2, 1, figsize=dim, dpi=400, sharex=True,sharey=True)

                astroPloting.radiiThickThin(ax, ax1, xaxis, mean_radii_Thin, mean_radii_Thick,
                                            poi, conv_1_style, r_outer_style, flux_peak_style, action,
                                            blurr_policy=True)

                figname = radVNu_path + model + "_Radii.jpeg"
                plt.savefig(figname, bbox_inches='tight')
                print("Image '{}' Created".format(figname))
                plt.close()

                # --------------------------------------------------
                # --------------------------------------------------
                # --------------------------------------------------

                '''Full Images----------------------------------'''
                parent_model_path = self.sub_paths["intensityPath"] + model + "/"
                current_model_file = parent_model_path + "blurr" + str(blur_kernal) + "/"

                done_list = None
                if blurr_frequency_list is not None:
                    done_list = np.full(len(blurr_frequency_list), False)

                k = action["start"]
                L = 0
                print("Constructing Full images for " + model)
                for i in range(num_of_action_points):
                    brightparams = self.all_model_brightparams[j]
                    brightparams[action["var"]] = k

                    current_freqeuncy = brightparams[action["var"]]

                    do_image, done_list = fileloading.frequencyListAnalysis(blurr_frequency_list, done_list,
                                                                            current_freqeuncy)

                    if do_image:
                        print("Full image production for intensity frame: ", L)
                        print(R"Observation frequency $\nu=$", str(k/1e9) + " GHz")

                        blurr_intensity_path = (current_model_file +
                                                action["var"] + "_blurr_" + "{:.5e}".format(current_freqeuncy))

                        lim0 = 25

                        print("Reading file: ", blurr_intensity_path)

                        h5f = h5py.File(blurr_intensity_path, 'r')

                        thin_blurr_image = h5f['thin_blurr_image'][:]
                        Absorbtion_Image = h5f["thick_blurr_image"][:]

                        h5f.close()

                        thin_intensity = [thin_blurr_image]
                        thick_intensity = [Absorbtion_Image]
                        thin_radii = [radii_I0_Thin[L, :]]
                        thick_radii = [radii_FullAbsorption_Thick[L, :]]

                        fig, (ax0, ax1) = plt.subplots(1, 2, figsize=[15, 7], dpi=400)

                        astroPloting.fullImage(fig, ax0, ax1, lim0, thin_intensity, thick_intensity,
                                               thin_radii, thick_radii, theta, blurr_policy=True)

                        ax0.text(-9, 8.5, var_label[action["var"]]
                                 + str(round(x_variable[L] / scale_label[action["var"]], 2))
                                 + ' ' + units_label[action["var"]], fontsize=12, color="w")

                        pltname = (image_path + '_FullImage_' + str(L) + "_Nu_"
                                   + str(round(x_variable[L] / scale_label[action["var"]], 2)) + ".jpeg")
                        plt.savefig(pltname, bbox_inches='tight')
                        print("Jpeg Created:  " + pltname)
                        plt.close()


                        '''IntensityVSRadii________________________________________________________________'''
                        fig, dum = plt.subplots(1, 2, figsize=[15, 7], dpi=400)
                        ax0 = plt.subplot(1, 2, 1)
                        ax1 = plt.subplot(1, 2, 2)
                        rmax = Absorbtion_Image.shape[0] * .4
                        astroPloting.IntensityVSRadii(fig, ax0, ax1, params.limits, thin_intensity,
                                                      thick_intensity, rmax,blurr_policy=True)
                        #
                        # ax1.text(2, 1.01, var_label[action["var"]]
                        #          + str(round(x_variable[i] / scale_label[action["var"]], 2))
                        #          + ' ' + units_label[action["var"]], fontsize=12, color="k")

                        pltname = (fluxVRadii_path + 'IntVRad_' + str(i) + "_Nu_"
                                   + str(round(x_variable[L] / scale_label[action["var"]], 2)) + ".jpeg")
                        plt.savefig(pltname, bbox_inches='tight')
                        print("Image '{}' Created".format(pltname))
                        plt.close()

                        '''RadVSVarphiType2________________________________________________________________'''
                        fig, dum = plt.subplots(1, 2, figsize=[15, 7], dpi=400)
                        ax0 = plt.subplot(1, 2, 1)
                        ax1 = plt.subplot(1, 2, 2)

                        astroPloting.radiiVSVarphi(fig, ax0, ax1, params.limits, thin_intensity,average=True,blurr_policy=True)

                        # ax0.text(2, 1.01, var_label[action["var"]]
                        #          + str(round(x_variable[i] / scale_label[action["var"]], 2))
                        #          + ' ' + units_label[action["var"]], fontsize=12, color="k")

                        pltname = (radVVarphi_path + 'radVVarphu_' + str(i) + "_Nu_"
                                   + str(round(x_variable[L] / scale_label[action["var"]], 2)) + ".jpeg")
                        plt.savefig(pltname, bbox_inches='tight')
                        print("Image '{}' Created".format(pltname))
                        plt.close()

                        L += 1
                    else:
                        print("Freuquency {} marked for skipping...".format("{:.5e}".format(current_freqeuncy)))
                    k += action['step']

            else:
                print(model + " marked for skipping...")

            j += 1  # marker for which brightparams to use

        if do_list is None:
            # histograms
            print(line)
            print("Creating Histograms")
            peak_hist_thin_path = self.sub_paths["peakHistThin"] + "blurr" + str(blur_kernal)
            peak_hist_thick_path = self.sub_paths["peakHistThick"] + "blurr" + str(blur_kernal)
            conv_hist_path = self.sub_paths["convHist"] + "blurr" + str(blur_kernal)
            total_flux_path = self.sub_paths["totalFlux"]+ "blurr" + str(blur_kernal)

            bar_xaxis = np.arange(len(self.all_model_names))
            """Flux Peaks_____________________________________________________________________"""
            # Thin_________________

            fig, ax = plt.subplots(1, 1, figsize=dim, dpi=400)

            astroPloting.histogram(ax, hist_flux_peaks_thins, "Flux Peak location (GHz)",
                                   "Optically Thin Assumption Frequency")

            figname = peak_hist_thin_path + "_FluxPeakThin.jpeg"
            plt.savefig(figname, bbox_inches='tight')
            print("Image '{}' Created".format(figname))
            plt.close()

            # Thick_______________

            fig, ax = plt.subplots(1, 1, figsize=dim, dpi=400)

            astroPloting.histogram(ax, hist_flux_peaks_thicks, "Flux Peak location (GHz)", "Full Solution Frequency")

            figname = peak_hist_thick_path + "_FluxPeakThick.jpeg"
            plt.savefig(figname, bbox_inches='tight')
            print("Image '{}' Created".format(figname))
            plt.close()

            # Thin Bar__________________
            fig, ax = plt.subplots(1, 1, figsize=dim, dpi=400)

            astroPloting.bar(ax, bar_xaxis, hist_flux_peaks_thins, "Optical Thin Assumption Peak Flux per model",
                             "Observation Frequency (GHz)", self.all_model_names)

            figname = peak_hist_thin_path + "_FluxPeakPerThinModel.jpeg"
            plt.savefig(figname, bbox_inches='tight')
            print("Image '{}' Created".format(figname))
            plt.close()

            # Thick Bar_______________

            fig, ax = plt.subplots(1, 1, figsize=dim, dpi=400)

            astroPloting.bar(ax, bar_xaxis, hist_flux_peaks_thicks, "Optical Thin Assumption Peak Flux per model",
                             "Observation Frequency (GHz)", self.all_model_names)

            figname = peak_hist_thick_path + "_FluxPeakPerThickModel.jpeg"
            plt.savefig(figname, bbox_inches='tight')
            print("Image '{}' Created".format(figname))
            plt.close()
            # Thick Bar

            """Conv_____________________________________________________________________"""
            # Conv Hist ______________________
            fig, ax = plt.subplots(1, 1, figsize=dim, dpi=400)

            astroPloting.histogram(ax, hist_convs, "Flux Peak location (GHz)", "Full Solution Frequency")

            figname = conv_hist_path + "_convHistFreqPerObservationFreq.jpeg"
            plt.savefig(figname, bbox_inches='tight')
            print("Image '{}' Created".format(figname))
            plt.close()

            # Conv Bar______________________
            fig, ax = plt.subplots(1, 1, figsize=dim, dpi=400)

            astroPloting.bar(ax, bar_xaxis, hist_convs, "Convergence for cumulative on I2",
                             "Observation Frequency (GHz)", self.all_model_names)
            figname = conv_hist_path + "_convHistPerModel.jpeg"
            plt.savefig(figname, bbox_inches='tight')
            print("Image '{}' Created".format(figname))
            plt.close()

            """230 total flux Thin___________________________________________________________"""
            fig, ax = plt.subplots(1, 1, figsize=dim, dpi=400)

            astroPloting.bar(ax, bar_xaxis, thin_total_flux, "Total Flux at 230GHz",
                             "Optically Thin Assumption Frequency", self.all_model_names)

            figname = total_flux_path + "_thin.jpeg"
            plt.savefig(figname, bbox_inches='tight')
            print("Image '{}' Created".format(figname))
            plt.close()

            """230 total flux Thick___________________________________________________________"""
            fig, ax = plt.subplots(1, 1, figsize=dim, dpi=400)

            astroPloting.bar(ax, bar_xaxis, thick_total_flux, "Total Flux at 230GHz",
                             "Full Solution Frequency", self.all_model_names)

            figname = total_flux_path + "_thick.jpeg"
            plt.savefig(figname, bbox_inches='tight')
            print("Image '{}' Created".format(figname))
            plt.close()
        else:
            print("Not all models were computed, skipping Histogram creation...")
    
    def getIntentPath(self,model,freq,blurr=None):
        if blurr is None:
            path = self.sub_paths['intensityPath'] + model + "/clean/nu0_" + "{:.5e}".format(freq)
        else:
            path = self.sub_paths['intensityPath'] + model + "/" + blurr +  "nu0_"+ "{:.5e}".format(freq)
        return path
    
    def getLBandPath(self,model):
        geo_name,dummy = self.getGeoName(model)
        return self.sub_paths['GeoDoth5Path'] + geo_name + "Lensing.h5"
    
    def getRBandPath(self,model):
        geo_name,dummy = self.getGeoName(model)
        return self.sub_paths['GeoDoth5Path'] + geo_name + "RayTracing.h5"
    
    def getPlotName(self,model):
        hlist = self.getHumanNameList(model)
        title = ''
        latexList = {
            'p_temp':   R"$\alpha_T$",
            'p_dens':   R"$\alpha_n$",
            'p_mag':    R"$\alpha_B$",
            'a':        R'$a$',
            'n_th0':    R'$n_{th,0}$',
            'theta_b':  R'$\theta_b$',
            'nu0':      R'$\nu_0$',
            'nscale':   R"$n_{scale}$",
            't_e0':     R'$t_{e,0}$',
            'b_0':      R'$b_0$'
        }
        for i in range(len(hlist)):
            if i == 0:
                pstr = R"${} = {}, "
            elif i == len(hlist) - 1:
                pstr= R"{} = {} $"
            else:   
                pstr = R"{} = {}, "
            title += pstr.format(latexList[hlist[i][0]],hlist[i][1])
        
        return title

    def getAllModelDict(self):
        diction = {}
        for i in range(len(self.all_model_names)):
            current_name = self.all_model_names[i]
            diction[current_name] = self.createModelDictionary(current_name)
        
        return diction
    
    def createModelDictionary(self,model):
        current_name = model
        geo_name,dummy = self.getGeoName(current_name)

        diction = {}
        diction['model'] = model
        diction['lband'] = self.getLBandPath(current_name)
        diction['rband'] = self.getRBandPath(current_name)
        diction['loadGeo'] = partial(fileloading.loadGeoModel,geo_name, self.run)
        diction['cintent_data'] = partial(self.getIntentPath,current_name)
        diction['pintent_data'] = self.sub_paths['intensityPath'] + current_name
        diction['bp'] = self.getModelBrightParams(current_name)
        diction['hname'] = self.makeHumanName(current_name)
        diction['hname_list'] =  self.getHumanNameList(current_name)
        diction['hname_plot'] = self.getPlotName(current_name)
        diction['np'] = {
            "jansky_thick": partial(self.getModelNp,current_name,'janksys_thick.npy'),
            "jansky_thin": partial(self.getModelNp,current_name,'janksys_thin.npy'),
            "radVnu_thick": partial(self.getModelNp,current_name,'mean_radii_Thick.npy'),
            "radVnu_thin": partial(self.getModelNp,current_name,'mean_radii_Thin.npy'),
            "xvar": partial(self.getModelNp,current_name,'x_variable.npy'),
            "jansky_thick": partial(self.getModelNp,current_name,'janksys_thick.npy'),
            "tau0": partial(self.getModelNp,current_name,'mean_optical_depth_I0.npy'),
            "tau1": partial(self.getModelNp,current_name,'mean_optical_depth_I1.npy'),
            "tau2": partial(self.getModelNp,current_name,'mean_optical_depth_I2.npy'),
            "tauT": partial(self.getModelNp,current_name,'mean_optical_depth_Total.npy')
        }
        diction['rprofs'] = partial(self.getRadialProfiles,model,None,None,(2,20),'log',500)
        diction['flux_peak'] = partial(self.getFluxPeak,current_name)
        diction['conv'] = partial(self.getConv,model)
        diction['image'] = partial(self.getIntensityName,model)

        t = os.listdir(diction['pintent_data'] + '/clean')
        t.remove('numpy')
        for k in range(len(t)):
            t[k] = float(t[k].replace("nu0_",""))
        t.sort()
        diction['cfreq_list'] = t

        return diction
    
    def getIntensityName(self, model,frequency,blurred=False,kernal=None):
        parent_model_path = self.sub_paths["intensityPath"] + model + "/"
        if not blurred:
            blurr_str = "clean/nu0_"
        else:
            if kernal is None:
                raise ValueError("No blurring kernal specified")
            blurr_str = "blurr" + str(kernal) + "/nu0_blurr_"

        return parent_model_path + blurr_str + "{:.5e}".format(frequency) 
    
    def getConv(self,model,method=1):
        yvar1 = self.getModelNp(model,'mean_radii_Thick.npy')[:,2]
        yvar2 = self.getModelNp(model,'mean_radii_Thick.npy')[:,3]
        xvar = self.getModelNp(model,'x_variable.npy')

        methods = {
            1: ilp.ring_convergance(xvar,yvar1,yvar2,2),
            2: ilp.ring_diff_min(xvar,yvar1,yvar2),
            3: ilp.ring_conv_3(xvar,yvar1,yvar2,2)
        }
        return methods[method]
    
    def getFluxPeak(self,model,ring,blurr=False,kernal=None):
        yvar = self.getModelNp(model,'janksys_thick.npy',blurr,kernal)[:,ring]
        xvar = self.getModelNp(model,'x_variable.npy',blurr,kernal)
        return self.getxOfFuncPeak(xvar,yvar)

    def getxOfFuncPeak(self,xvar,yvar):
        num_of_observation_points = 1000
        
        coords = np.linspace(xvar.min(),xvar.max(),num_of_observation_points)
        points = num_of_observation_points - 1
        step = (xvar.max() - xvar.min())/points
        
        yvar_interp = interp1d(xvar,yvar)(coords)
        return (xvar.min() + step * np.argmax(yvar_interp))
        

    def getModelNp(self,model,array:str,blurr=False,kernal=None):
        if blurr is True:
            if kernal is None:
                raise ValueError("No blurring kernal specified")
            numpy_str = "/blurr" + str(kernal) + "/numpy/blurr_"
        else:
            numpy_str = "/clean/numpy/"
        path = self.sub_paths['intensityPath'] + model + numpy_str + array
        return np.load(path)
    
    def getModelGrid(self):
        dimensions = [0]
        for i in range(len(self.geo_grid_names)):
            dimensions[0] += 1

        for i in range(self.run_type):
            current_dimension = self.var_intensity_grid_names[i]
            dimensions += [self.variable_param_ranges[current_dimension]]
        
        return [self.oneGridDepth,self.twoGridDepth,self.threeGridDepth][self.run_type - 1](dimensions)
        
    def oneGridDepth(self,dimensions):
        grid = np.ndarray([dimensions[0],dimensions[1]],dtype=object)
        l = 0
        for i in range(dimensions[0]):
            for j in range(dimensions[1]):
                grid[i,j] = self[str(self.all_model_names[l])]
                l +=1
        
        return grid
    
    def twoGridDepth(self,dimensions):
        grid = np.ndarray([dimensions[0],dimensions[1],dimensions[2]],dtype=object)
        l = 0
        for i in range(dimensions[0]):
            for j in range(dimensions[1]):
                for k in range(dimensions[2]):
                    grid[i,j,k] = self[str(self.all_model_names[l])]
                    l +=1
        return grid
    
    def threeGridDepth(self,dimensions):
        grid = np.ndarray([dimensions[0],dimensions[1],dimensions[2],dimensions[3]],dtype=object)
        l = 0
        for i in range(dimensions[0]):
            for j in range(dimensions[1]):
                for k in range(dimensions[2]):
                    for p in range(dimensions[3]):
                        grid[i,j,k,p] = self[str(self.all_model_names[l])]
                        l +=1
        return grid


    def __getitem__(self,item):
        if type(item) is str:
            return self.createModelDictionary(item)
        elif type(item) is tuple:
            return self.getModelGrid()[item]

    def getGeoName(self,model):
        modelstr = model.split("_")
        if len(modelstr) == 1:
            modelstr = modelstr[0]
            splitchar = ''
            l = 0
            while not splitchar.isnumeric():
                splitchar = modelstr[l]
                l += 1
            modelstr = model.split(splitchar)
            geo_name = modelstr[0]
            int_name = model.split(geo_name[len(geo_name)-1])[1]
        else:
            geo_name = modelstr[0]
            int_name = modelstr[1]
        return geo_name, int_name

    def makeHumanName(self,model):
        
        geo_name, int_name = self.getGeoName(model)
        geo_param = self.geo_grid_params[self.geo_grid_names.index(geo_name)]
        # pstr = 
        hname = ''
        for j in range(len(geo_param[0])):
            hname += '--' + geo_param[0][j] + ">" + geo_param[1][j]
        hname = hname.removeprefix("--")
            
        for j in range(len(self.var_intensity_grid_names)):
            current_var_param = self.var_intensity_grid_names[j]
            current_var_position = int(int_name[j]) - 1

            value = str(self.string_order[current_var_param].format(self.intensity_grid_params[current_var_param][current_var_position]))
            hname += "--{}>{}".format(current_var_param,value)
        
        return hname
    
    def getHumanNameList(self,model):
        geo_name, int_name = self.getGeoName(model)
        geo_param = self.geo_grid_params[self.geo_grid_names.index(geo_name)]
        # pstr = 
        hname = []
        for j in range(len(geo_param[0])):
            hname += [(geo_param[0][j], geo_param[1][j])]
        
            
        for j in range(len(self.var_intensity_grid_names)):
            current_var_param = self.var_intensity_grid_names[j]
            current_var_position = int(int_name[j]) - 1

            value = str(self.string_order[current_var_param].format(self.intensity_grid_params[current_var_param][current_var_position]))
            hname += [(current_var_param,value)]
        return hname

    def getRadialProfiles(self,model,rmin=None,rmax=None,rrange=None,scale:str='lin',num_of_points=100,theta_bkey=funckeys["theta_bkey"],nu0=230e9,lab_frame=False):
        if ((rmin is not None) and (rmax is None)) or ((rmax is not None) and (rmin is None)):
            raise ValueError("Either rmin or rmax were specified while the other was")
        if rrange is None:
            if (rmax is None) or (rmin is None):
                raise ValueError("either rmin and rmax or rrange must be specified")
        
        if (rrange is not None) and (rmin is not None):
            raise ValueError("rrange and rmin/rmax can't be simulaneously specified")
        if (scale != 'lin') and (scale != 'log'):
            raise ValueError("Scale '{}' unrecongnized".format(scale))
        
        if rmin is None:
            rmin = rrange[0]
            rmax = rrange[1]
        
        # _______________
        bp = self.getModelBrightParams(model)
        bp = bp.copy()
        bp["nu0"]=nu0
        bp = self.add_units_to_bp(bp)
        if scale == 'lin':
            rarray = np.linspace(rmin,rmax,num_of_points)
        elif scale == 'log':
            rarray = np.logspace(np.log10(rmin),np.log10(rmax),num_of_points)

        b_profile = ilp.b_func_power_variable(rarray,mass=bp['mass'],rb_0=bp['rb_0'],bv_0=bp['b_0'],p_bv=bp['p_mag'])
        te_profile = ilp.te_func(rarray,mass=bp['mass'],rb_0=bp['rb_0'],t_e0=bp['t_e0'], p_temp=bp['p_temp'])

        theta_e=ilp.theta_e_func(te_profile)
        nth_profile = ilp.nth_func(rarray,mass=bp['mass'],rb_0=bp['rb_0'], n_th0=bp['n_th0'], p_dens=bp['p_dens'])
        
        if theta_bkey is None:
            theta_bkey = self.funckeys["theta_bkey"]

        funckeys_use = self.funckeys.copy()
        funckeys_use["theta_bkey"] = theta_bkey
        si_thin, si_thick, d, d = ilp.thermal_profile({'r':rarray,'x':0,'y':0}, redshift=1, cosAng=1,bp=bp, fk=funckeys_use)

        profiles = {
            'r':rarray,
            'b': b_profile,
            'theta_e':theta_e,
            'nth':nth_profile,
            'si_thick': si_thick,
            'si_thin': si_thin
                    }
        
        if lab_frame is True:
            final_data_path = self.sub_paths["intensityPath"] + model + "/" + "clean/" + "numpy/"
            
            profiles['n0_si_thick'] = np.load(final_data_path + "_full_profiles0_{}GHz".format("{:.5e}".format(nu0)) + ".npy")
            profiles['n1_si_thick'] = np.load(final_data_path + "_full_profiles1_{}GHz".format("{:.5e}".format(nu0)) + ".npy")
            profiles['n2_si_thick'] = np.load(final_data_path + "_full_profiles2_{}GHz".format("{:.5e}".format(nu0)) + ".npy")
        return profiles



    def add_units_to_bp(self,bp:dict):
        bp_unit = bp.copy()
        bp_unit['nu0'] = bp['nu0'] * ilp.Hz
        bp_unit['mass'] = bp['mass'] * u.g
        bp_unit['theta_b'] = bp['theta_b'] * ilp.rads
        bp_unit['n_th0'] = bp['n_th0'] * ilp.cmcubed
        bp_unit['t_e0'] = bp['t_e0'] * ilp.kelv
        bp_unit['b_0'] = bp['b_0'] * ilp.gauss
        return bp_unit

        

def fmt(x, pos):
    x = x / 1e9
    return '{:.2f}'.format(x)
