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
import re

breaker = "     "
line = "\n________________________\n"
long_line = "\n________________________________________________________________________\n"
line_small = "________________________________ \n"


class BigRuns:

    def __init__(self, run: str, intensity_grid_params: dict, var_intensity_grid_names, geo_grid_params, geo_grid_names,
                 normalized_brightparams=False):
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
            "3d": image_path + "3dPloting/"

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

        print("Final Grid Type: ",self.run_type)
        # Create astro ParamGrid________________________________________________________________________________________
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
                arrays[k] = np.load(file,allow_pickle=True)
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
                2: self.type2Grid
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

    def type2Grid(self):
        all_brightparams = []  # list[tuple (name, bright parameters)]
        all_model_names = []
        for i in range(self.variable_param_ranges[self.var_intensity_grid_names[0]]):
            for j in range(self.variable_param_ranges[self.var_intensity_grid_names[1]]):
                current_model_brightparams = {}  # brightprams for current model
                current_model_name = "Model" + str(i + 1) + str(j + 1)
                # Fill out the constant parameters
                for key in list(self.constant_params):
                    current_model_brightparams[key] = self.constant_params[key]

                current_model_brightparams[self.var_intensity_grid_names[0]] =\
                    self.intensity_grid_params[self.var_intensity_grid_names[0]][i]
                current_model_brightparams[self.var_intensity_grid_names[1]] =\
                    self.intensity_grid_params[self.var_intensity_grid_names[1]][j]

                all_brightparams += [current_model_brightparams]
                all_model_names += [current_model_name]

        return all_model_names,all_brightparams

    def type1Grid(self):
        all_brightparams= []
        all_model_names = []
        for i in range(self.variable_param_ranges[self.var_intensity_grid_names[0]]):
            current_model_brightparams = {}  # dict[key] = for value of [key]/varied parameter
            current_model_name = "Model" + str(i + 1)
            for key in list(self.constant_params):
                current_model_brightparams[key] = self.constant_params[key]

            current_model_brightparams[self.var_intensity_grid_names[0]] = self.intensity_grid_params[self.var_intensity_grid_names[0]][i]

            all_brightparams += [current_model_brightparams]
            all_model_names += [current_model_name]

        return all_model_names,all_brightparams

    def type0grid(self):
        current_model_brightparams = {}  # dict[key] = for value of [key]/varied parameter
        current_model_name = "Model" + str(1)
        for key in list(self.constant_params):
            current_model_brightparams[key] = self.constant_params[key]

        return [current_model_name], [current_model_brightparams]

    def createGeoGrid(self):
        for i in range(len(self.geo_grid_names)):
            # Create normalizing lensing band

            normal_param_name = self.geo_grid_names[i] + "_Normalizing"

            fileloading.loadGeoModel(normal_param_name, self.run)
            importlib.reload(params)

            new_lband = self.sub_paths["GeoDoth5Path"] + normal_param_name + "Lensing" + ".h5"
            new_rtray = self.sub_paths["GeoDoth5Path"] + normal_param_name + "RayTracing" + ".h5"
            print("Computation of {} Lensing Bands".format(normal_param_name))
            """Computation of the lensing bands______________________________________________________"""
            subprocess.run(['python3 ' + EZPaths.aartPath + '/lensingbands.py '], shell=True)

            spin_case = params.spin_case
            i_case = params.i_case

            fnbands = path + "LensingBands_a_%s_i_%s.h5" % (params.spin_case, params.i_case)

            print("Reading file: ", fnbands)

            '''Analytical Ray-tracing______________________________________________________'''
            print("Analytical Raytracing of {}".format(normal_param_name))
            subprocess.run(['python3 ' + EZPaths.aartPath + '/raytracing.py '], shell=True)
            fnrays1 = path + "Rays_a_%s_i_%s.h5" % (spin_case, i_case)

            # Move lensing bands and raytracing bands
            subprocess.run(["mv " + fnbands + ' ' + new_lband], shell=True)
            subprocess.run(["mv " + fnrays1 + ' ' + new_rtray], shell=True)

            # Actual Model Params____________________________________
            fileloading.loadGeoModel(self.geo_grid_names[i], self.run)
            importlib.reload(params)

            new_lband = self.sub_paths["GeoDoth5Path"] + self.geo_grid_names[i] + "Lensing" + ".h5"
            new_rtray = self.sub_paths["GeoDoth5Path"] + self.geo_grid_names[i] + "RayTracing" + ".h5"
            print("Computation of {} Lensing Bands".format(self.geo_grid_names[i]))
            """Computation of the lensing bands______________________________________________________"""
            subprocess.run(['python3 ' + EZPaths.aartPath + '/lensingbands.py '], shell=True)

            spin_case = params.spin_case
            i_case = params.i_case

            fnbands = path + "LensingBands_a_%s_i_%s.h5" % (params.spin_case, params.i_case)

            print("Reading file: ", fnbands)

            '''Analytical Ray-tracing______________________________________________________'''
            print("Analytical Raytracing of {}".format(self.geo_grid_names[i]))
            subprocess.run(['python3 ' + EZPaths.aartPath + '/raytracing.py '], shell=True)
            fnrays1 = path + "Rays_a_%s_i_%s.h5" % (spin_case, i_case)

            # Move lensing bands and raytracing bands
            subprocess.run(["mv " + fnbands + ' ' + new_lband], shell=True)
            subprocess.run(["mv " + fnrays1 + ' ' + new_rtray], shell=True)

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
                self.all_model_brightparams += [self.all_intensity_model_brightparams[i]]
                self.total_models_count += 1
                k += 1

    def createBrightParamsDict(self):
        diction = {}
        for i in range(len(self.all_model_brightparams)):
            # String Names
            diction[self.all_model_names[i]] = self.all_model_brightparams[i]
            # ________________________________
        return diction

    def creatIntensityGrid(self,action,do_list=None,isContinuous=False):

        funckeys = {
            "emodelkey": 0,  # emodelkey Emission Model choice, 0 = thermal ultrarelativistic, 1 = power law
            "bkey": 2,  # bkey
            "nnoisykey": 0,  # nnoisykey Inoisy density. 0 = no noise, 1 = noise
            "tnoisykey": 0,  # tnoisykey Inoisy temperature
            "bnoisykey": 0  # bnoisykey Inoisy magnetic field
        }

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

            for i in range(len(self.all_intensity_model_brightparams)):
                # String Names
                current_total_name = self.all_model_names[k]
                # ________________________________
                current_bp = self.all_model_brightparams[k]

                current_model_file = self.sub_paths["intensityPath"] + current_total_name + "/clean/"

                # only skips below if run is continuous and file already exist, or if do is false
                preform_model = fileloading.crossContinousDoAnalysis(
                    current_total_name,do_list,current_model_file,isContinuous)
                if preform_model:
                    if not self.already_normalized_brightparams:
                        print("\n" + "Normalizing " + current_total_name + "\n")
                        print(long_line)

                        current_bp["n_th0"] = normalizingBrightparams.normalize(normlband, normrtray, current_bp)

                        print("\n" + current_total_name + " normalized with a value of n_th0="
                              + str(current_bp["n_th0"]) + "\n")

                    print("Creating Intensity Movie for Model ", current_total_name)
                    print(long_line)

                    if self.run_type == 0:
                        run_type_arg = 1
                    else:
                        run_type_arg = self.run_type
                    intermodel_data = movieMakerIntensity.intensity_movie(
                        action, self.sub_paths, current_total_name, run_type_arg, current_bp)

                    print(
                        "\nTotal flux at 230GHz for Optically Thin Assumption: " + str(intermodel_data["thin_total_flux"]))
                    print("Total flux at 230GHz for Full Solution: " + str(intermodel_data["thick_total_flux"]) + "\n")
                    all_230_total_jy_thin += [intermodel_data["thin_total_flux"]]
                    all_230_total_jy_thick += [intermodel_data["thick_total_flux"]]
                else:
                    print("{} marked for skipping...".format(current_total_name))
                k += 1
        if not self.already_normalized_brightparams:
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

        self.already_normalized_brightparams = True
        # Numpy saving________________________________________________

        all_230_total_jy_thin_numpy_name = self.sub_paths["meta"] + "thin_total_flux"
        all_230_total_jy_thick_numpy_name = self.sub_paths["meta"] + "thick_total_flux"

        np.save(all_230_total_jy_thin_numpy_name, np.array(all_230_total_jy_thin))
        np.save(all_230_total_jy_thick_numpy_name, np.array(all_230_total_jy_thick))

    def intensityGridAnalysis(self,action,do_list=None,isContinuous=False):
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
                print("Analyzing Intensity Movie for Model ", current_total_name)
                print(long_line)
                movieMakerIntensity.imageAnalysis(
                    action, self.sub_paths, current_total_name, current_bp
                )
            else:
                print(current_total_name + " marked for skipping...")

    def blurrIntensityGrid(self,action,do_list=None,isContinuous=False):

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
            current_model_file = parent_model_path + "blurr/"

            preform_model = fileloading.crossContinousDoAnalysis(
                current_total_name, do_list, current_model_file, isContinuous)

            if preform_model:
                intermodel_data = movieMakerIntensity.blurr_intensity_movie(
                    action, self.sub_paths, current_total_name, run_type_arg, current_bp)

                print(
                    "\nTotal flux at 230GHz for Optically Thin Assumption: " + str(intermodel_data["thin_total_flux"]))
                print("Total flux at 230GHz for Full Solution: " + str(intermodel_data["thick_total_flux"]) + "\n")
                all_230_total_jy_thin += [intermodel_data["thin_total_flux"]]
                all_230_total_jy_thick += [intermodel_data["thick_total_flux"]]
            else:
                print(current_total_name + " marked for skipping...")
        # Numpy saving________________________________________________

        all_230_total_jy_thin_numpy_name = self.sub_paths["meta"] + "blurr_thin_total_flux"
        all_230_total_jy_thick_numpy_name = self.sub_paths["meta"] + "blurr_thick_total_flux"

        np.save(all_230_total_jy_thin_numpy_name, np.array(all_230_total_jy_thin))
        np.save(all_230_total_jy_thick_numpy_name, np.array(all_230_total_jy_thick))

    def blurrIntensityGridAnalysis(self,action,do_list=None,isContinuous=False):
        print(line)
        print(line)
        print(line)
        print("Analyzing blurred intensity grid for " + self.run)
        for i in range(len(self.all_model_brightparams)):

            # String Names
            current_total_name = self.all_model_names[i]
            # ________________________________
            current_bp = self.all_model_brightparams[i]

            parent_model_path = self.sub_paths["intensityPath"] + current_total_name + "/"
            current_model_file = parent_model_path + "blurr/"

            preform_model = fileloading.crossContinousDoAnalysis(
                current_total_name, do_list, current_model_file, isContinuous)

            if preform_model:
                print("Analyzing blurred intensity movie for model ", current_total_name)

                movieMakerIntensity.blurrImageAnalysis(
                    action, self.sub_paths, current_total_name, current_bp
                )
            else:
                print(current_total_name + " marked for skipping...")

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
                intensity_model_string += breaker + key + ": " + str(self.constant_params[key]) + '\n'

        intensity_model_string += line_small
        for i in range(len(self.all_model_names)):
            current_name = self.all_model_names[i]
            current_model = self.all_model_brightparams[i]
            intensity_model_string += line_small + current_name + '\n'

            for k in range(len(self.var_intensity_grid_names)):
                intensity_model_string += (breaker + self.var_intensity_grid_names[k] + ": "
                                           + str(current_model[self.var_intensity_grid_names[k]]) + '\n')
            intensity_model_string += breaker + "n_th0: " + str(current_model["n_th0"]) + '\n'

        intensity_model_string += line_small + line_small + line_small

        geo_models_string = "\nGEOMODELS\n"
        for i in range(len(self.geo_grid_names)):
            geo_models_string += self.geo_grid_names[i] + "| "
            for j in range(len(self.geo_grid_params[i][0])):
                geo_models_string += '\n    ' + self.geo_grid_params[i][0][j] + ": " + self.geo_grid_params[i][1][j]
            geo_models_string += '\n'

        return intensity_model_string + geo_models_string

    def graphCreation(self,action,do_list=None,isContinuous=False):
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
        hist_flux_peaks_thins = []
        hist_flux_peaks_thicks = []
        hist_convs = []
        dim = [10, 8]
        for model in self.all_model_names:

            fluxVNu_path = self.sub_paths["fluxPath"] + model + "/clean/"
            radVNu_path = self.sub_paths["radPath"] + model + "/clean/"
            image_path = self.sub_paths["imagePath"] + model + "/clean/"
            Optical_depth_path = self.sub_paths["opticalDepth"] + model + "/clean/"

            radVVarphi_path = self.sub_paths["RadVVarphi"] + model + "/clean/"
            fluxVRadii_path = self.sub_paths["fluxVRadii"] + model + "/clean/"

            file_creation = [fluxVNu_path, radVNu_path, image_path, Optical_depth_path, radVVarphi_path,
                             fluxVRadii_path]

            preform_model = fileloading.crossContinousDoAnalysis(
                model, do_list, fluxVNu_path, isContinuous)

            if preform_model:
                print(line)
                print("Graphing " + model)
                for i in range(len(file_creation)):
                    fileloading.creatSubDirectory(file_creation[i], kill_policy=True)

                fluxVNu_path += "Clean_"
                radVNu_path += "Clean_"
                image_path += "Clean_"
                Optical_depth_path += "Clean_"
                radVVarphi_path += "Clean_"
                fluxVRadii_path += "Clean_"


                if self.run_type == 0:
                    amount_to_subtract = 1
                else:
                    amount_to_subtract = self.run_type

                current_geo_model = model[0:len(model) - amount_to_subtract]
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

                # Points of Interest

                conv_1 = (action["start"] + action["step"] *
                          ilp.ring_convergance(mean_radii_Thick[:, 2], mean_radii_Thick[:,3],3))
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
                if action["images"]:
                    parent_model_path = self.sub_paths["intensityPath"] + model + "/"
                    current_model_file = parent_model_path + "clean/"
                    k = action["start"]
                    print("Constructing Full images for " + model)
                    for i in range(num_of_intensity_points):

                        brightparams = self.all_model_brightparams[j]
                        brightparams["nu0"] = k
                        print("Full image production for intensity frame: ", i)
                        print(R"Observation frequency $\nu=$",k)

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
                        thin_intensity = [I0,I1,I2,I0 + I1 + I2]
                        thick_intensity = [I0_Absorb, I1_Absorb,I2_Absorb, Absorbtion_Image]
                        thin_radii = [radii_I0_Thin[i,:],radii_I1_Thin[i,:],radii_I2_Thin[i,:],radii_Full_Thin[i,:]]
                        thick_radii = [radii_I0_Thick[i, :], radii_I1_Thick[i, :],
                                       radii_I2_Thick[i, :], radii_FullAbsorption_Thick[i, :]]

                        fig, (ax0, ax1) = plt.subplots(1, 2, figsize=[15, 7], dpi=400)

                        astroPloting.fullImage(fig,ax0,ax1,lim0,thin_intensity, thick_intensity,
                                               thin_radii, thick_radii, theta)

                        ax0.text(-9, 8.5, astroModels.var_label[action["var"]]
                                 + str(round(x_variable[i] / astroModels.scale_label[action["var"]], 2))
                                 + ' ' + astroModels.units_label[action["var"]], fontsize=12, color="w")

                        pltname = (image_path + 'FullImage_' + str(i) + "_Nu_"
                                   + str(round(x_variable[i] / astroModels.scale_label[action["var"]], 2)) + ".jpeg")
                        plt.savefig(pltname, bbox_inches='tight')
                        print("Jpeg Created:  " + pltname)
                        plt.close()

                        # Get total jansky
                        # VARPHI__________________________________________________-

                        '''IntensityVSRadiiType1________________________________________________________________'''
                        fig, dum = plt.subplots(1, 2, figsize=[15, 7], dpi=400)
                        ax0 = plt.subplot(1, 2, 1)
                        ax1 = plt.subplot(1, 2, 2)

                        astroPloting.IntensityVSRadiiType1(fig, ax0, ax1, params.limits, thin_intensity,
                                                           thick_intensity, rmax)
                        #
                        # ax1.text(2, 1.01, astroModels.var_label[action["var"]]
                        #          + str(round(x_variable[i] / astroModels.scale_label[action["var"]], 2))
                        #          + ' ' + astroModels.units_label[action["var"]], fontsize=12, color="k")

                        pltname = (fluxVRadii_path + 'IntVRad_' + str(i) + "_Nu_"
                                   + str(round(k / astroModels.scale_label[action["var"]], 2)) + ".jpeg")
                        plt.savefig(pltname, bbox_inches='tight')
                        print("Image '{}' Created".format(pltname))
                        plt.close()

                        '''RadVSVarphiType2________________________________________________________________'''
                        fig, dum = plt.subplots(1, 2, figsize=[15, 7], dpi=400)
                        ax0 = plt.subplot(1, 2, 1)
                        ax1 = plt.subplot(1, 2, 2)

                        astroPloting.radiiVSVarphi(fig, ax0, ax1, params.limits, thin_intensity)

                        # ax0.text(2, 1.01, astroModels.var_label[action["var"]]
                        #          + str(round(x_variable[i] / astroModels.scale_label[action["var"]], 2))
                        #          + ' ' + astroModels.units_label[action["var"]], fontsize=12, color="k")

                        pltname = (radVVarphi_path + 'radVVarphu_' + str(i) + "_Nu_"
                                   + str(round(k / astroModels.scale_label[action["var"]], 2)) + ".jpeg")
                        plt.savefig(pltname, bbox_inches='tight')
                        print("Image '{}' Created".format(pltname))
                        plt.close()


                        k += action['step']
                j += 1  # marker for which brightparams to use
            else:
                print(model + " marked for skipping...")
        # histograms
        print(line)
        print("Creating Histograms")
        peak_hist_thin_path = self.sub_paths["peakHistThin"]
        peak_hist_thick_path = self.sub_paths["peakHistThick"]
        conv_hist_path = self.sub_paths["convHist"]
        total_flux_path = self.sub_paths["totalFlux"]

        bar_xaxis = np.arange(len(self.all_model_names))
        bar_labels = self.all_model_names
        for i in range(len(bar_xaxis)):
            bar_labels[i] = bar_labels[i].replace("Model","")
        """Flux Peaks_____________________________________________________________________"""
        # Thin_________________

        fig, ax = plt.subplots(1, 1, figsize=dim, dpi=400)

        astroPloting.histogram(ax,hist_flux_peaks_thins,"Flux Peak location (GHz)","Optically Thin Assumption Frequency")

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

        astroPloting.bar(ax,bar_xaxis,hist_flux_peaks_thins,"Optical Thin Assumption Peak Flux per model",
                         "Observation Frequency (GHz)",bar_labels)

        figname = peak_hist_thin_path + "FluxPeakPerThinModel.jpeg"
        plt.savefig(figname, bbox_inches='tight')
        print("Image '{}' Created".format(figname))
        plt.close()

        # Thick Bar_______________

        fig, ax = plt.subplots(1, 1, figsize=dim, dpi=400)

        astroPloting.bar(ax,bar_xaxis,hist_flux_peaks_thicks,"Optical Thin Assumption Peak Flux per model",
                         "Observation Frequency (GHz)",bar_labels)

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

        astroPloting.bar(ax,bar_xaxis,hist_convs,"Convergence for cumulative on I2",
                         "Observation Frequency (GHz)",bar_labels)
        figname = conv_hist_path + "convHistPerModel.jpeg"
        plt.savefig(figname, bbox_inches='tight')
        print("Image '{}' Created".format(figname))
        plt.close()

        """230 total flux Thin___________________________________________________________"""
        fig, ax = plt.subplots(1, 1, figsize=dim, dpi=400)

        astroPloting.bar(ax, bar_xaxis,thin_total_flux, "Total Flux at 230GHz",
                         "Optically Thin Assumption Frequency",bar_labels)

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


    def blurrGraphCreation(self,action):
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

            if self.run_type == 0:
                amount_to_subtract = 1
            else:
                amount_to_subtract = self.run_type

            current_geo_model = model[0:len(model) - amount_to_subtract]
            fileloading.loadGeoModel(current_geo_model, self.run)

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
            data_path = self.sub_paths["intensityPath"] + model + "/blurr/numpy/"

            x_variable = np.load(data_path + "blurr_x_variable.npy")
            janksys_thick = np.load(data_path + "blurr_janksys_thick.npy")
            janksys_thin = np.load(data_path + "blurr_janksys_thin.npy")
            mean_radii_Thin = np.load(data_path + "blurr_mean_radii_Thin.npy")
            mean_radii_Thick = np.load(data_path + "blurr_mean_radii_Thick.npy")
            radii_I0_Thin = np.load(data_path + "blurr_radii_I0_Thin.npy")
            radii_FullAbsorption_Thick = np.load(data_path + "blurr_radii_FullAbsorption_Thick.npy")
            theta = np.load(data_path + "blurr_theta.npy")

            num_of_intensity_points = janksys_thin[:,0].shape[0]
            print("Number of Intensity Points: ", num_of_intensity_points)

            xaxis = np.array(x_variable) / astroModels.scale_label[action['var']]
            # one_M = ilp.rg_func(brightparams["mass"] * u.g).to(u.m)
            # M2uas = np.arctan(one_M.value / dBH) / muas_to_rad

            fluxVNu_path = self.sub_paths["fluxPath"] + model + "/blurr/"
            radVNu_path = self.sub_paths["radPath"] + model + "/blurr/"
            image_path = self.sub_paths["imagePath"] + model + "/blurr/"
            Optical_depth_path = self.sub_paths["opticalDepth"] + model + "/blurr/"

            file_creation = [fluxVNu_path, radVNu_path,image_path,Optical_depth_path]

            for i in range(len(file_creation)):
                fileloading.creatSubDirectory(file_creation[i],kill_policy=True)

            fluxVNu_path += "Blurr_"
            radVNu_path += "Blurr_"
            image_path += "Blurr_"
            Optical_depth_path += "Blurr_"

            # Points of Interest

            conv_1 = action["start"] + action["step"] * ilp.ring_convergance(mean_radii_Thick[:, 0], mean_radii_Thick[:, 0],
                                                                             3)
            conv_1 = conv_1 / astroModels.scale_label[action['var']]

            flux_peak_thin = action["start"] + action["step"] * np.argmax(janksys_thin[:, 0])
            flux_peak_thin = flux_peak_thin / astroModels.scale_label[action['var']]

            flux_peak_thick = action["start"] + action["step"] * np.argmax(janksys_thick[:, 0])
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

            # _______________________________________________
            # _______________________________________________
            # ________________________________
            '''JANKSKY PLOTS----------------------------------'''

            fig, (ax, ax1) = plt.subplots(2, 1, figsize=dim, dpi=400, sharex=True)

            astroPloting.fluxThickThin(ax, ax1, xaxis, janksys_thin, janksys_thick,
                                       poi, conv_1_style, r_outer_style, flux_peak_style, action,blurr_policy=True)

            figname = fluxVNu_path + model + "_Flux.jpg"
            plt.savefig(figname, bbox_inches='tight')
            plt.close()
            print("Image '{}' Created".format(figname))

            # ______________________________________________
            # ______________________________________________
            # ___________________________________
            '''RADII PLOTS----------------------------------'''
            fig, (ax, ax1) = plt.subplots(2, 1, figsize=dim, dpi=400, sharex=True)

            astroPloting.radiiThickThin(ax, ax1, xaxis, mean_radii_Thin, mean_radii_Thick,
                                        poi, conv_1_style, r_outer_style, flux_peak_style, action,blurr_policy=True)

            figname = radVNu_path + model + "_Radii.jpeg"
            plt.savefig(figname, bbox_inches='tight')
            print("Image '{}' Created".format(figname))
            plt.close()

            # --------------------------------------------------
            # --------------------------------------------------
            # --------------------------------------------------

            '''Full Images----------------------------------'''
            if action["images"]:
                parent_model_path = self.sub_paths["intensityPath"] + model + "/"
                current_model_file = parent_model_path + "blurr/"

                k = action["start"]
                print("Constructing Full images for " + model)
                for i in range(num_of_intensity_points):
                    brightparams = self.all_model_brightparams[j]
                    brightparams["nu0"] = k
                    print("Full image production for intensity frame: ", i)
                    print(R"Observation frequency $\nu=$",k)

                    blurr_intensity_path = (current_model_file +
                                            action["var"] + "_blurr_" + "{:.5e}".format(brightparams[action["var"]]))

                    lim0 = 25

                    print("Reading file: ", blurr_intensity_path)

                    h5f = h5py.File(blurr_intensity_path, 'r')

                    thin_blurr_image = h5f['thin_blurr_image'][:]
                    Absorbtion_Image = h5f["thick_blurr_image"][:]

                    h5f.close()

                    thin_intensity = [thin_blurr_image]
                    thick_intensity = [Absorbtion_Image]
                    print("SHAPE: ",radii_I0_Thin.shape)
                    thin_radii = [radii_I0_Thin[i,:]]
                    thick_radii = [radii_FullAbsorption_Thick[i, :]]

                    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=[15, 7], dpi=400)

                    astroPloting.fullImage(fig,ax0,ax1,lim0,thin_intensity, thick_intensity,
                                           thin_radii, thick_radii, theta,blurr_policy=True)

                    ax0.text(-9, 8.5, astroModels.var_label[action["var"]]
                             + str(round(x_variable[i] / astroModels.scale_label[action["var"]], 2))
                             + ' ' + astroModels.units_label[action["var"]], fontsize=12, color="w")

                    pltname = (image_path + '_FullImage_' + str(i) + "_Nu_"
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
        peak_hist_thin_path = self.sub_paths["peakHistThin"]
        peak_hist_thick_path = self.sub_paths["peakHistThick"]
        conv_hist_path = self.sub_paths["convHist"]
        total_flux_path = self.sub_paths["totalFlux"]

        bar_xaxis = np.arange(len(self.all_model_names))
        """Flux Peaks_____________________________________________________________________"""
        # Thin_________________

        fig, ax = plt.subplots(1, 1, figsize=dim, dpi=400)

        astroPloting.histogram(ax,hist_flux_peaks_thins,"Flux Peak location (GHz)","Optically Thin Assumption Frequency")

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

        astroPloting.bar(ax,bar_xaxis,hist_flux_peaks_thins,"Optical Thin Assumption Peak Flux per model",
                         "Observation Frequency (GHz)",self.all_model_names)

        figname = peak_hist_thin_path + "_FluxPeakPerThinModel.jpeg"
        plt.savefig(figname, bbox_inches='tight')
        print("Image '{}' Created".format(figname))
        plt.close()

        # Thick Bar_______________

        fig, ax = plt.subplots(1, 1, figsize=dim, dpi=400)

        astroPloting.bar(ax,bar_xaxis,hist_flux_peaks_thicks,"Optical Thin Assumption Peak Flux per model",
                         "Observation Frequency (GHz)",self.all_model_names)

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

        astroPloting.bar(ax,bar_xaxis,hist_convs,"Convergence for cumulative on I2",
                         "Observation Frequency (GHz)",self.all_model_names)
        figname = conv_hist_path + "_convHistPerModel.jpeg"
        plt.savefig(figname, bbox_inches='tight')
        print("Image '{}' Created".format(figname))
        plt.close()

        """230 total flux Thin___________________________________________________________"""
        fig, ax = plt.subplots(1, 1, figsize=dim, dpi=400)

        astroPloting.bar(ax, bar_xaxis,thin_total_flux, "Total Flux at 230GHz",
                         "Optically Thin Assumption Frequency",self.all_model_names)

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

    def smallModelRun(self,models:list[str],action,save_path):
        for model in models:
            print("Running " + model)
            print(line)
            print(line)
            print("Initializing Graph Creation")
            thin_total_flux = np.load(self.sub_paths["meta"] + "thin_total_flux.npy")
            thick_total_flux = np.load(self.sub_paths["meta"] + "thick_total_flux.npy")

            current_geo_model = fileloading.totalModelNametoGridModel(model,self.run_type)
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
            data_path = self.sub_paths["intensityPath"] + model + "/" + "numpy/"

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

            num_of_intensity_points = int((action["stop"] - action["start"]) / action["step"])
            print("Number of Intensity Points: ", num_of_intensity_points)

            print("Constructing Full images for " + model)
            if action["images"]:
                k = action["start"]
                for i in range(num_of_intensity_points):
                    print(self.all_model_names)
                    brightparams = self.all_model_brightparams[self.all_model_names.index(model)]
                    brightparams["nu0"] = k
                    print("Full image production for intensity frame: ", i)
                    print(R"Observation frequency $\nu=$", k)

                    current_intensity_file = (self.sub_paths["intensityPath"] + model + "/" + action["var"]
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

                    astroPloting.IntensityVSRadiiType1(fig, ax0, ax1, ax2, ax3, params.limits, thin_intensity,
                                                       thick_intensity, rmax)
                    #
                    # ax1.text(2, 1.01, astroModels.var_label[action["var"]]
                    #          + str(round(x_variable[i] / astroModels.scale_label[action["var"]], 2))
                    #          + ' ' + astroModels.units_label[action["var"]], fontsize=12, color="k")

                    pltname = (save_path['intVRad'] + 'IntVRad_' + str(i) + "_Nu_"
                               + str(round(k / astroModels.scale_label[action["var"]], 2)) + ".jpeg")
                    plt.savefig(pltname, bbox_inches='tight')
                    print("Image '{}' Created".format(pltname))
                    plt.close()

                    '''IntensityVSRadiiType2________________________________________________________________'''
                    fig, dum = plt.subplots(1, 2, figsize=dim, dpi=400)
                    ax0 = plt.subplot(1, 2, 1)
                    ax1 = plt.subplot(1, 2, 2)

                    astroPloting.IntensityVSRadiiType2(fig, ax0, ax1, params.limits, thin_intensity, rmax)

                    # ax0.text(2, 1.01, astroModels.var_label[action["var"]]
                    #          + str(round(x_variable[i] / astroModels.scale_label[action["var"]], 2))
                    #          + ' ' + astroModels.units_label[action["var"]], fontsize=12, color="k")

                    pltname = (save_path['intVRad2'] + 'IntVRad_' + str(i) + "_Nu_"
                               + str(round(k / astroModels.scale_label[action["var"]], 2)) + ".jpeg")
                    plt.savefig(pltname, bbox_inches='tight')
                    print("Image '{}' Created".format(pltname))
                    plt.close()

                    '''RadVSVarphiType2________________________________________________________________'''
                    fig, dum = plt.subplots(1, 2, figsize=dim, dpi=400)
                    ax0 = plt.subplot(1, 2, 1)
                    ax1 = plt.subplot(1, 2, 2)

                    astroPloting.radiiVSVarphi(fig, ax0, ax1, params.limits, thin_intensity)

                    # ax0.text(2, 1.01, astroModels.var_label[action["var"]]
                    #          + str(round(x_variable[i] / astroModels.scale_label[action["var"]], 2))
                    #          + ' ' + astroModels.units_label[action["var"]], fontsize=12, color="k")

                    pltname = (save_path['radVVarphi'] + 'radVVarphu_' + str(i) + "_Nu_"
                               + str(round(k / astroModels.scale_label[action["var"]], 2)) + ".jpeg")
                    plt.savefig(pltname, bbox_inches='tight')
                    print("Image '{}' Created".format(pltname))
                    plt.close()


                    '''RadVSVarphiType2Intensities________________________________________________________________'''
                    # fig, dum = plt.subplots(1, 2, figsize=dim, dpi=400)
                    # ax0 = plt.subplot(1, 2, 1)
                    # ax1 = plt.subplot(1, 2, 2)
                    #
                    # astroPloting.radiiVSVarphi(fig, ax0, ax1, params.limits, thin_intensity,plot_intensites=True)
                    #
                    # # ax0.text(2, 1.01, astroModels.var_label[action["var"]]
                    # #          + str(round(x_variable[i] / astroModels.scale_label[action["var"]], 2))
                    # #          + ' ' + astroModels.units_label[action["var"]], fontsize=12, color="k")
                    #
                    # pltname = (save_path['radVVarphi'] + 'radVVarphuIntesities_' + str(i) + "_Nu_"
                    #            + str(round(k / astroModels.scale_label[action["var"]], 2)) + ".jpeg")
                    # plt.savefig(pltname, bbox_inches='tight')
                    # print("Image '{}' Created".format(pltname))
                    # plt.close()

                    k += action['step']


def fmt(x, pos):
    x = x / 1e9
    return '{:.2f}'.format(x)


