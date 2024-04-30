import sys
import subprocess

aartpath = '/home/td6241/repositories/aart' #insert path to aart repo
sys.path.append(aartpath)

import shutil
import subprocess
import EZPaths
import os
from aart_func import *
import params
from params import *
import importlib
import astroModels


def createIntensityArgs(brightparams,funckeys=astroModels.funckeys):
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


def intensityNameWrite(brightparams,funckeys):
    filename = path + ('Intensity_a_{}_i_{}_nu_{}_mass_{}_scaleh_{}_thetab_{}_beta_{}_rie_{}_rb_{}_nth0_{}_te0_{}_'
                       'b0_{}_pdens_{}_ptemp_{}_pmag_{}_nscale_{}_emkey_{}_bkey_{}_nkey_{}_tnkey_{}_bnkey_{}.h5').format(
        params.spin_case,
        i_case,
        "{:.5e}".format(brightparams["nu0"].value),
        "{:.5e}".format(brightparams["mass"].value),
        float(brightparams["scale_height"]),
        "{:.3f}".format(brightparams["theta_b"].value),
        "{:.2f}".format(float(brightparams["beta"])),
        "{:.1f}".format(float(brightparams["r_ie"])),
        "{:.1f}".format(float(brightparams["rb_0"])),
        "{:.1e}".format(brightparams["n_th0"].value),
        "{:.1e}".format(brightparams["t_e0"].value),
        "{:.3e}".format(brightparams["b_0"].value),
        float(brightparams["p_dens"]),
        float(brightparams["p_temp"]),
        float(brightparams["p_mag"]),
        "{:.1f}".format(brightparams["nscale"]),
        funckeys["emodelkey"],
        funckeys["bkey"],
        funckeys["nnoisykey"],
        funckeys["tnoisykey"],
        funckeys["bnoisykey"]
    )

    return filename


def intensityNameNoUnits(brightparams,funckeys):

    filename = path + ('Intensity_a_{}_i_{}_nu_{}_mass_{}_scaleh_{}_thetab_{}_beta_{}_rie_{}_rb_{}_nth0_{}_te0_{}_'
                       'b0_{}_pdens_{}_ptemp_{}_pmag_{}_nscale_{}_emkey_{}_bkey_{}_nkey_{}_tnkey_{}_bnkey_{}.h5').format(
        params.spin_case,
        i_case,
        "{:.5e}".format(brightparams["nu0"]),
        "{:.5e}".format(brightparams["mass"]),
        float(brightparams["scale_height"]),
        "{:.3f}".format(brightparams["theta_b"]),
        "{:.2f}".format(float(brightparams["beta"])),
        "{:.1f}".format(float(brightparams["r_ie"])),
        "{:.1f}".format(float(brightparams["rb_0"])),
        "{:.1e}".format(brightparams["n_th0"]),
        "{:.1e}".format(brightparams["t_e0"]),
        "{:.3e}".format(brightparams["b_0"]),
        float(brightparams["p_dens"]),
        float(brightparams["p_temp"]),
        float(brightparams["p_mag"]),
        "{:.1f}".format(brightparams["nscale"]),
        funckeys["emodelkey"],
        funckeys["bkey"],
        funckeys["nnoisykey"],
        funckeys["tnoisykey"],
        funckeys["bnoisykey"]
    )

    return filename


def loadGeoModel(current_model:str, run:str):
    print("\nLoading " + current_model + "\n")
    paramFile = EZPaths.aartPath + "/params.py"
    model_is_loaded = os.path.exists(paramFile)
    if model_is_loaded:
        subprocess.run(["rm " + paramFile], shell=True)

    loading_model = EZPaths.modelRunsDir + run + "/geoModels/" + current_model + ".py"
    loaded_model = EZPaths.aartPath + "/" + "params.py"
    cmd = "cp " + loading_model + " " + loaded_model
    subprocess.run([cmd], shell=True)
    print("file {} Loaded as {}".format(loading_model,loaded_model))

    importlib.reload(params)


def numOfModel(model:str,all_full_model_names):
    j = 0
    for i in range(len(all_full_model_names)):
        if all_full_model_names[i] == model:
            print(all_full_model_names[i] + " Identified as " + model + " for number " + str(j))
            return j
        j += 1


def creatSubDirectory(path_chosen,purpose:str = '',kill_policy=False):
    isDir = os.path.exists(path_chosen)
    if kill_policy:
        if isDir:
            print("Subdirectory for " + purpose + " '{}' already exist, removing...".format(path_chosen))

            existing_file = os.listdir(path_chosen)
            if len(existing_file) == 0:
                print("Empty directory, removing...")
                os.rmdir(path_chosen)
            else:
                print("Not empty directory, removing...")
                shutil.rmtree(path_chosen)

            os.makedirs(path_chosen)
            print("A Subdirectory for " + purpose + " '{}' was created".format(path_chosen))
        else:
            os.makedirs(path_chosen)
            print("A Subdirectory for " + purpose + " '{}' was created".format(path_chosen))

    else:
        if not isDir:
            os.makedirs(path_chosen)
            print("A Subdirectory for " + purpose + " '{}' was created".format(path_chosen))
        else:
            print("Subdirectory for " + purpose + " '{}' already exist, doing nothing".format(path_chosen))


def runsInit(run:str, grid_params:dict, var_params):
    """
    Args:
        run: name of this run
        grid_params: dictionary containing all intensity parameters. Each value is a
                     list containing all parameters to be computed in current run
        var_params: list of key names to be varable parameters in the current run,
                    correspond to more than 1 entry in list for grid params
    Returns:

    """

    creatSubDirectory(EZPaths.modelRunsDir, "storing all run results")
    main_path = EZPaths.modelRunsDir + run + "/"
    creatSubDirectory(main_path, "storing this run result")

    image_path = main_path + "Images/"
    data_path = main_path + "Data/"
    meta_path = main_path + "InterModel/"
    creatSubDirectory(image_path,"images")
    creatSubDirectory(data_path, "data")
    creatSubDirectory(meta_path, "inter model data")

    sub_paths = {
        "GeoDoth5Path": data_path + "geo/",
        "intensityPath": data_path + "intensity/",
        "fluxPath": image_path + "fluxVNu/",
        "radPath":image_path + "radVNu/",
        "imagePath":image_path + "image/",
        "opticalDepth": image_path + "opticalDepth/",
        "peakHistThin": image_path + "peakHistThin/",
        "peakHistThick": image_path + "peakHistThick/",
        "convHist": image_path + "convHist/",
        "totalFlux":image_path + "totalFlux230Ghz/",
        "meta": meta_path,
        "3d": image_path + "3dPloting/"

    }

    for key in list(sub_paths):
        # creatSubDirectory(sub_paths[key])
        isDir = os.path.exists(sub_paths[key])
        if not isDir:
            os.makedirs(sub_paths[key])
            print("Subdirectory {} Created".format(sub_paths[key]))

    all_intensity_models,total_models_count, run_type, variable_param_ranges, constant_param= (
        createAstroParams(grid_params, var_params))

    return sub_paths, all_intensity_models, total_models_count, run_type, variable_param_ranges, constant_param


def createAstroParams(intensity_grid_params:dict, var_params:list[str]):
    """

    Args:
        intensity_grid_params: dictionary containing all intensity parameters. Each value is a
                     list containing all parameters to be computed in current run
        var_params: list of key names to be varable parameters in the current run,
                    correspond to more than 1 entry in list for grid params

    Returns: list[tuple(model_name,brightparams)] without Normalization

    """

    total_models_count, run_type, variable_param_ranges, constant_param = (
        analyzeRunType(intensity_grid_params, var_params)
    )

    # Create the Models
    grid_types = {
        2: type2Grid
    }

    all_models = grid_types[run_type](intensity_grid_params, var_params,
                                      variable_param_ranges, constant_param, total_models_count)

    return all_models, total_models_count, run_type, variable_param_ranges, constant_param


def analyzeRunType(intensity_grid_params:dict, var_params:list[str]):
    """

    Args:
        intensity_grid_params: dictionary containing all intensity parameters. Each value is a
                     list containing all parameters to be computed in current run
        var_params: list of key names to be varable parameters in the current run,
                    correspond to more than 1 entry in list for grid params

    Returns:

    """

    total_models_count = 1
    run_type = 0
    variable_param_ranges = {}
    constant_param = {}

    for key in list(intensity_grid_params):

        param_range = len(intensity_grid_params[key])
        # print("Key, " + key + " count: " + str(var_params.count(key)))
        if param_range > 1 and var_params.count(key) == 0:
            raise ValueError("More than one parameter put in for static parameter")

        if var_params.count(key) == 1:
            run_type += 1
            variable_param_ranges[key] = param_range
            total_models_count = total_models_count * param_range
        else:
            constant_param[key] = intensity_grid_params[key][0]

    return total_models_count, run_type, variable_param_ranges, constant_param


def type2Grid(grid_params:dict, var_params:list[str], variable_param_ranges:dict,
              constant_params:dict, total_models_count:int):
    """

    Args:
        grid_params: dictionary containing all intensity parameters. Each value is a
                     list containing all parameters to be computed in current run
        var_params: list of key names to be varable parameters in the current run,
                      correspond to more than 1 entry in list for grid params
        variable_param_ranges: contain number of values to run through for a given key/parameter
        constant_params: parameters that will remain constant for all models
        total_models_count: total number of models for this run

    Returns: list[tuple(model_name,brightparams)]

    """
    all_models = []  # list[tuple (name, bright parameters)]
    for i in range(variable_param_ranges[var_params[0]]):
        for j in range(variable_param_ranges[var_params[1]]):
            current_model = {}  # brightparams
            current_model_name = "Model" + str(i + 1) + str(j + 1)
            # Fill out the constant parameters
            for key in list(constant_params):
                current_model[key] = constant_params[key]

            current_model[var_params[0]] = grid_params[var_params[0]][i]
            current_model[var_params[1]] = grid_params[var_params[1]][j]

            all_models += [(current_model_name, current_model)]

    return all_models


line = "________________________________ \n"
breaker = "     "


def tupleToString(all_models:list[tuple], var_params, constant_params, total_models_count):
    """

    Args:
        all_models: list[(name:str, bright parameters:dict)]
        var_params: list of key names to be varable parameters in the current run,
              correspond to more than 1 entry in list for grid params
        constant_params: parameters that will remain constant for all models
        total_models_count: total number of models for this run

    Returns:

    """
    string = line + line + line + "Total Number of Models: " + str(total_models_count) + '\n' + "Constant Params: " + '\n'
    for key in list(constant_params):
        string += breaker + key + ": " + str(constant_params[key]) + '\n'

    string += line
    for i in range(len(all_models)):
        string += line + all_models[i][0] + '\n'
        current_model = all_models[i][1]
        for k in range(len(var_params)):
            string += breaker + var_params[k] + ": " + str(current_model[var_params[k]]) + '\n'

    return string + line + line + line


def totalModelNametoGridModel(model,run_type):
    if run_type == 0:
        subtract_amount = 1
    else:
        subtract_amount = run_type
    return model[0:len(model) - subtract_amount]


def writeDocString(filename:str, text:str):
    doc_string_file = filename
    cmd = "touch " + doc_string_file
    print("File ({}) Created".format(filename))
    subprocess.run([cmd], shell=True)

    # writing
    doc_string_file = open(doc_string_file, 'w')
    # doc_string_file.write(full_string + geo_models_string)
    doc_string_file.write(text)
    doc_string_file.close()


def doListAnalysis(current_model_name, do_list):
    do = False
    if current_model_name in do_list:
        do = True

    return do


def continousAnalysis(file_name,isContinuous):
    # returns False if and only if is continous and the file already exist
    return not (isContinuous and os.path.exists(file_name))


def crossContinousDoAnalysis(current_model_name, do_list, file_name, isContinuous):
    # If "do_list" is not Null, only factor in continuous condition and "do" returns false no matter what
    # If "do_list" is none, only factor in continous condition
    if do_list is not None:
        return doListAnalysis(current_model_name, do_list)
    else:
        return continousAnalysis(file_name,isContinuous)


def frequencyListAnalysis(frequency_list, done_list, current_frequency):
    """

    Args:
        frequency_list: list of intensity points to blurr at
        done_list: list of bool values indicating which blurr points have already been done
        current_frequency: current frequency

    Returns: updated blurr list, extra True value if condition met

    """
    first_false = np.argmin(done_list)
    do_image = False

    if (((frequency_list is not None) and (current_frequency >= frequency_list[first_false]))
            and (not done_list[first_false])):  # if all done, min gives first argument
        done_list[first_false] = True
        do_image = True

    if frequency_list is None:
        do_image = True

    return do_image,done_list







