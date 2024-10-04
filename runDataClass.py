import astroModels
import numpy as np
import movieMakerV2.movieMakerIntensity
from movieMakerV2 import movieMakerIntensity
from astroModels import *

class runData:

    def __init__(self, name: str, brightparams_grid: dict, variable_brightparams_names: list,
                 variable_geo_values, variable_geo_names, funckey=funckeys):
        """

        Args:
            name: Name of run
            brightparams_grid: dict[parameter A] = [A1,A2,...AN]
            variable_brightparams_names: list of all variables whose list in the brightparam grid is longer than 1
            variable_geo_values: list[ModelA_list, ModelB_list,...] where
                ModelA_tuple([paramA,paramB,...], [ParamA_value,paramB_value,...])
            variable_geo_names:  list of geometrical file name
        """
        self.run_name = name
        self.brightparams_values_grid = brightparams_grid
        self.variable_brightparams_names = variable_brightparams_names
        self.variable_geo_values_grid = variable_geo_values
        self.variable_geo_names = variable_geo_names
        self.action = {
            "var": "nu0",
            "start": 10e9,
            "stop": 700e9,
            "step": 20e9,
            "images": True
        }

        self.isNormalized = False
        self.blurr_policy = False
        self.funckeys = funckey

    def setAction(self, action):
        """

        Args:
            action: {
                "var": str,
                "start": Float,
                "stop": Float,
                "step": Float,
                "images": Bolean

        """
        self.action = action

    def setActionKey(self, key, value):
        self.action[key] = value

    def setisNormalized(self, isnormalized):
        self.isNormalized = isnormalized

    def setBlurrPolicy(self, blurr_policy):
        self.blurr_policy = blurr_policy

    def getRunName(self):
        return self.run_name

    def getBrightparams(self):
        return self.brightparams_values_grid

    def getFunckeys(self):
        return self.funckeys

    def getBPVarNames(self):
        return self.variable_brightparams_names

    def getGeoGrid(self):
        return self.variable_geo_values_grid

    def getGeoGridNames(self):
        return self.variable_geo_names

    def getAction(self):
        return self.action

    def getIsNormalized(self):
        return self.isNormalized

    def getBlurrPolicy(self):
        return self.blurr_policy


# run1______________________________________________________________________
run1 = runData("run1",
               astroModels.bp_run1,
               ["p_temp", "p_mag"],
               [(["a"], [str(.3)]), (["a"], [str(.9)])],
               ["ModelA", "ModelB"])

# run2______________________________________________________________________
run2 = runData("run2",
               astroModels.bp_run2,
               ["p_temp", "p_mag"],
               [(["a"], [str(.001)]), (["a"], [str(.5)]), (["a"], [str(15 / 16)])],
               ["ModelA", "ModelB", "ModelC"],
               )
run2.setisNormalized(True)

# TestRun__________________________________________
TestRun  = runData("TestRun",
               astroModels.bp_PrunSingle,
               ["t_e0"],
               [(["a"], [str(15 / 16)])],
               ["ModelC"],
               funckey=astroModels.funckey_PrunSingle
               )
TestRun.setisNormalized(False)


# PrunSingle______________________________________________________________________
PrunSingle  = runData("PrunSingle",
               astroModels.bp_PrunSingle,
               [],
               [(["a"], [str(15 / 16)])],
               ["ModelC"],
               funckey=astroModels.funckey_PrunSingle
               )
PrunSingle.setisNormalized(True)


# PRUN______________________________________________________________________
PRUN  = runData("PRUN",
               astroModels.bp_PRUN,
               ["p_temp", "p_mag"],
               [(["a"], [str(.001)]), (["a"], [str(.5)]), (["a"], [str(15 / 16)])],
               ["ModelA", "ModelB", "ModelC"],
               )
PRUN.setisNormalized(True)

# exp1__________________________________________________________________________
exp1 = runData("exp1",
               astroModels.bp_exp1,
               ["t_e0"],
               [(["a"], [str(15 / 16)])],
               ["ModelC"],
               )
exp1.setisNormalized(True)

# exp2__________________________________________________________________________
exp2 = runData("exp2",
               astroModels.bp_exp2,
               ["t_e0", "b_0"],
               [(["a"], [str(15 / 16)])],
               ["ModelA"],
               )
exp2.setisNormalized(False)
exp2Action = {
    "var": "nu0",
    "start": 10e9,
    "stop": 100e10,
    "step": 40e9,
    "images": True
}
exp2.setAction(exp2Action)

# testRun1______________________________________________________________________

testRun1 = runData("testRun1",
                   astroModels.bp_testRun1,
                   ["p_temp", "p_mag"],
                   [(["a"], [str(.3)]), (["a"], [str(.9)])],
                   ["ModelA", "ModelB"])

testRun1.setAction(
    {
        "var": "nu0",
        "start": 3.00e+10,
        "stop": 6.00e+10,
        "step": 1.00e+10,
        "images": True
    }
)
# testRun2______________________________________________________________________
testRun2 = runData("testRun2",
                   astroModels.bp_testRun1,
                   ["p_temp", "p_mag"],
                   [(["a"], [str(.3)]), (["a"], [str(.9)])],
                   ["ModelA", "ModelB"])

testRun2.setAction(
    {
        "var": "nu0",
        "start": 3.00e+10,
        "stop": 6.00e+10,
        "step": 1.00e+10,
        "images": True
    }
)
testRun2.setisNormalized(True)

# soloRun1______________________________________________________________________
soloRun1 = runData("soloRun1",
                   astroModels.bp_soloRun1,
                   ["p_mag"],
                   [(["a"], [str(.3)]), (["a"], [str(.9)])],
                   ["ModelA", "ModelB"])

# soloRun2______________________________________________________________________
soloRun2 = runData("soloRun2",
                   astroModels.bp_soloRun2,
                   [],
                   [(["a"], [str(15 / 16)])],
                   ["ModelB"])
soloRun2.setAction(
    {
        "var": "nu0",
        "start": 670e9,
        "stop": 700e9,
        "step": 20e9,
        "images": True
    }
)
soloRun2.setisNormalized(True)


class SingleModelData:

    def __init__(self, sub_paths: dict, run: str, model: str, frequency_points=None):
        freq_points = frequency_points
        self.sub_paths = sub_paths
        self.run = run
        self.model_name = model
        self.average = True

        data_path = self.sub_paths["intensityPath"] + model + "/clean/numpy/"
        self.clean_data_paths = {
            "x_variable": data_path + "x_variable.npy",
            "janksys_thick": data_path + "janksys_thick.npy",
            "janksys_thin": data_path + "janksys_thin.npy",
            "mean_radii_Thin": data_path + "mean_radii_Thin.npy",
            "mean_radii_Thick": data_path + "mean_radii_Thick.npy",
            "radii_I0_Thin": data_path + "radii_I0_Thin.npy",
            "radii_I1_Thin": data_path + "radii_I1_Thin.npy",
            "radii_I2_Thin": data_path + "radii_I2_Thin.npy",
            "radii_Full_Thin": data_path + "radii_Full_Thin.npy",
            "radii_FullAbsorption_Thick": data_path + "radii_FullAbsorption_Thick.npy",
            "radii_I0_Thick": data_path + "radii_I0_Thick.npy",
            "radii_I1_Thick": data_path + "radii_I1_Thick.npy",
            "radii_I2_Thick": data_path + "radii_I2_Thick.npy",
            "theta": data_path + "theta.npy",
            "mean_optical_depth_I0": data_path + "mean_optical_depth_I0.npy",
            "mean_optical_depth_I1": data_path + "mean_optical_depth_I1.npy",
            "mean_optical_depth_I2": data_path + "mean_optical_depth_I2.npy",
            "mean_optical_depth_Total": data_path + "mean_optical_depth_Total.npy",
        }
        if freq_points is not None:
            for i in range(len(freq_points)):
                keys = [str(i) + "full_profiles0", str(i) + "full_profiles1", str(i) + "full_profiles2"]
                L = 0
                for key in keys:
                    self.clean_data_paths[key] = data_path + ("_full_profiles" + str(L) +
                                                              "_{}GHz.npy").format("{:.5e}".format(freq_points[i]))
                    L += 1

        data_path = self.sub_paths["intensityPath"] + model + "/blurr/numpy/"

        self.blurr_data_paths = {
            "x_variable": data_path + "blurr_x_variable.npy",
            "janksys_thick": data_path + "blurr_janksys_thick.npy",
            "janksys_thin": data_path + "blurr_janksys_thin.npy",
            "mean_radii_Thin": data_path + "blurr_mean_radii_Thin.npy",
            "mean_radii_Thick": data_path + "blurr_mean_radii_Thick.npy",
            "radii_I0_Thin": data_path + "blurr_radii_I0_Thin.npy",
            "radii_FullAbsorption_Thick": data_path + "blurr_radii_FullAbsorption_Thick.npy",
            "theta": data_path + "theta.npy",
            # "86full_profiles": data_path + "_full_profiles0_{}GHz.npy".format("{:.5e}".format(freq_points[0])),
            # # "86full_profile_units": data_path + "_full_profiles_unit_{}GHz.npy".format("{:.5e}".format(freq_points[0])),
            # "230full_profiles": data_path + "_full_profiles0_{}GHz.npy".format("{:.5e}".format(freq_points[1])),
            # # "230full_profile_units": data_path + "_full_profiles_unit_{}GHz.npy".format("{:.5e}".format(freq_points[1])),
            # "345full_profiles": data_path + "_full_profiles0_{}GHz.npy".format("{:.5e}".format(freq_points[2])),
            # # "345full_profile_units": data_path + "_full_profiles_unit_{}GHz.npy".format("{:.5e}".format(freq_points[2]))
        }

        self.give_clean = True

    def set_give_policy(self, policy: bool):
        self.give_clean = policy

    def set_average_policy(self, truth):
        self.average = truth
        data_path = self.sub_paths["intensityPath"] + self.model_name + "/clean/numpy/"

        if not self.average:
            data_path += "FalseAvg_"

        self.clean_data_paths = {
            "x_variable": data_path + "x_variable.npy",
            "janksys_thick": data_path + "janksys_thick.npy",
            "janksys_thin": data_path + "janksys_thin.npy",
            "mean_radii_Thin": data_path + "mean_radii_Thin.npy",
            "mean_radii_Thick": data_path + "mean_radii_Thick.npy",
            "radii_I0_Thin": data_path + "radii_I0_Thin.npy",
            "radii_I1_Thin": data_path + "radii_I1_Thin.npy",
            "radii_I2_Thin": data_path + "radii_I2_Thin.npy",
            "radii_Full_Thin": data_path + "radii_Full_Thin.npy",
            "radii_FullAbsorption_Thick": data_path + "radii_FullAbsorption_Thick.npy",
            "radii_I0_Thick": data_path + "radii_I0_Thick.npy",
            "radii_I1_Thick": data_path + "radii_I1_Thick.npy",
            "radii_I2_Thick": data_path + "radii_I2_Thick.npy",
            "theta": data_path + "theta.npy",
            "mean_optical_depth_I0": data_path + "mean_optical_depth_I0.npy",
            "mean_optical_depth_I1": data_path + "mean_optical_depth_I1.npy",
            "mean_optical_depth_I2": data_path + "mean_optical_depth_I2.npy",
            "mean_optical_depth_Total": data_path + "mean_optical_depth_Total.npy"
        }

    def __getitem__(self, item):
        if self.give_clean:
            return np.load(self.clean_data_paths[item], allow_pickle=True)
        else:
            return np.load(self.blurr_data_paths[item], allow_pickle=True)

    def get_intensity_name(self, frequency):
        parent_model_path = self.sub_paths["intensityPath"] + self.model_name + "/"
        current_model_file = parent_model_path + "clean/"
        return current_model_file + "nu0" + "_" + "{:.5e}".format(frequency)

    def get_sub_paths(self):
        return self.sub_paths
