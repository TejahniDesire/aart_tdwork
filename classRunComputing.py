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

line = "________________________________ \n"
breaker = "     "


class BigRuns:

    def __init__(self, run: str, grid_params: dict, var_params):
        """
        Args:
            run: name of this run
            grid_params: dictionary containing all intensity parameters. Each value is a
                         list containing all parameters to be computed in current run
            var_params: list of key names to be varable parameters in the current run,
                        correspond to more than 1 entry in list for grid params
        Returns:

        """
        self.intensity_grid_params = grid_params
        self.run = run
        self.var_params = var_params

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
            "fluxPath": image_path + "fluxVNu/",
            "radPath": image_path + "radVNu/",
            "imagePath": image_path + "image/",
            "opticalDepth": image_path + "opticalDepth/",
            "peakHistThin": image_path + "peakHistThin/",
            "peakHistThick": image_path + "peakHistThick/",
            "convHist": image_path + "convHist/",
            "totalFlux": image_path + "totalFlux230Ghz/",
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
            if param_range > 1 and self.var_params.count(key) == 0:
                raise ValueError("More than one parameter put in for static parameter")

            if self.var_params.count(key) == 1:
                self.run_type += 1
                self.variable_param_ranges[key] = param_range
                self.total_intensity_models_count = self.total_intensity_models_count * param_range
            else:
                self.constant_params[key] = self.intensity_grid_params[key][0]

        # Create astro ParamGrid________________________________________________________________________________________
        grid_types = {
            0: self.type0grid,
            1: self.type1Grid,
            2: self.type2Grid
        }

        self.all_models = grid_types[self.run_type]()

        print(self.tupleToString())

    def type2Grid(self):
        all_models = []  # list[tuple (name, bright parameters)]
        for i in range(self.variable_param_ranges[self.var_params[0]]):
            for j in range(self.variable_param_ranges[self.var_params[1]]):
                current_model_brightparams = {}  # brightprams for current model
                current_model_name = "Model" + str(i + 1) + str(j + 1)
                # Fill out the constant parameters
                for key in list(self.constant_params):
                    current_model_brightparams[key] = self.constant_params[key]

                current_model_brightparams[self.var_params[0]] = self.intensity_grid_params[self.var_params[0]][i]
                current_model_brightparams[self.var_params[1]] = self.intensity_grid_params[self.var_params[1]][j]

                all_models += [(current_model_name, current_model_brightparams)]

        return all_models

    def type1Grid(self):
        all_models = []
        for i in range(self.variable_param_ranges[self.var_params[0]]):
            current_model_brightparams = {}  # dict[key] = for value of [key]/varied parameter
            current_model_name = "Model" + str(i + 1)
            for key in list(self.constant_params):
                current_model_brightparams[key] = self.constant_params[key]

            current_model_brightparams[self.var_params[0]] = self.intensity_grid_params[self.var_params[0]][i]

            all_models += [(current_model_name, current_model_brightparams)]

        return all_models

    def type0grid(self):
        current_model_brightparams = {}  # dict[key] = for value of [key]/varied parameter
        current_model_name = "Model" + str(1)
        for key in list(self.constant_params):
            current_model_brightparams[key] = self.constant_params[key]

        return [(current_model_name, current_model_brightparams)]

    def tupleToString(self):
        """

        Args:
            all_models: list[(name:str, bright parameters:dict)]
            var_params: list of key names to be varable parameters in the current run,
                  correspond to more than 1 entry in list for grid params
            constant_params: parameters that will remain constant for all models
            total_models_count: total number of models for this run

        Returns:

        """
        string = line + line + line + "Total Number of Models: " + str(
            self.total_intensity_models_count) + '\n' + "Constant Params: " + '\n'
        for key in list(self.constant_params):
            string += breaker + key + ": " + str(self.constant_params[key]) + '\n'

        string += line
        for i in range(len(self.all_models)):
            string += line + self.all_models[i][0] + '\n'
            current_model = self.all_models[i][1]
            for k in range(len(self.var_params)):
                string += breaker + self.var_params[k] + ": " + str(current_model[self.var_params[k]]) + '\n'

        return string + line + line + line
