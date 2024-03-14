import subprocess
import EZPaths
import os
from aart_func import *
import params
from params import *
import importlib
import astroModels


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
    print("Loading " + current_model)
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


def runsInit(run:str, grid_params, var_params):
    """
    var_params: list of key names to be varable parameters in the current run
    """
    isDir = os.path.exists(EZPaths.modelRunsDir)
    if not isDir:
        os.makedirs(EZPaths.modelRunsDir)
        print("A directory ({}) was created to store the run results".format(EZPaths.modelRunsDir))

    main_path = EZPaths.modelRunsDir + run + "/"
    isDir = os.path.exists(main_path)
    if not isDir:
        os.makedirs(main_path)
        print("A directory ({}) was created to store the run results".format(main_path))

    sub_paths = {
        "GeoDoth5Path": main_path + "geo/",
        "intensityPath": main_path + "intensity/",
        "fluxPath": main_path + "fluxVNu/",
        "radPath":main_path + "radVNu/",
        "imagePath":main_path + "image/",
        "opticalDepth": main_path + "opticalDepth/"
    }
    for key in list(sub_paths):
        isDir = os.path.exists(sub_paths[key])
        if not isDir:
            os.makedirs(sub_paths[key])
            print("Subdirectory {} Created".format(sub_paths[key]))

    all_models, names = createAstroParams(grid_params, var_params)

    return sub_paths, all_models, names


def createAstroParams(grid_params:dict,var_params:list[str]):

    # Read Inputted Models
    total_models_count = 1
    total_model_types = 0
    variable_param_ranges = {}
    constant_param = {}
    for key in list(grid_params):

        param_range = len(grid_params[key])
        # print("Key, " + key + " count: " + str(var_params.count(key)))
        if param_range > 1 and var_params.count(key) == 0:
            raise ValueError("More than one parameter put in for static parameter")

        if var_params.count(key) == 1:
            total_model_types += 1
            variable_param_ranges[key] = param_range
            total_models_count = total_models_count * param_range
        else:
            constant_param[key] = grid_params[key][0]

    # Create the Models
    grid_types = {
        2: type2Grid
    }

    all_models, names = grid_types[total_model_types](grid_params,var_params,variable_param_ranges, constant_param, total_models_count)

    return all_models, names


def type2Grid(grid_params:dict,defined_list, variable_param_ranges, constant_param, total_models_count):
    all_models = [] # tuple (name, bright parameters)
    for i in range(variable_param_ranges[defined_list[0]]):
        for j in range(variable_param_ranges[defined_list[1]]):
            current_model = {}
            current_model_name = "Model" + str(i + 1) + str(j + 1)
            # Fill out the constant parameters
            for key in list(constant_param):
                current_model[key] = constant_param[key]

            current_model[defined_list[0]] = grid_params[defined_list[0]][i]
            current_model[defined_list[1]] = grid_params[defined_list[1]][j]

            all_models += [(current_model_name, current_model)]

    return all_models, tupleToString(all_models, defined_list, constant_param,total_models_count)


line = "________________________________ \n"
breaker = "     "


def tupleToString(all_models, defined_list, constant_params, total_models_count):
    string = line + line + line + "Total Number of Models: " + str(total_models_count) + '\n' + "Constant Params: " + '\n'
    for key in list(constant_params):
        string += breaker + key + ": " + str(constant_params[key]) + '\n'

    string += line
    for i in range(len(all_models)):
        string += line + all_models[i][0] + '\n'
        current_model = all_models[i][1]
        for k in range(len(defined_list)):
            string += breaker + defined_list[k] + ": " + str(current_model[defined_list[k]]) + '\n'

    return string + line + line + line





# loadGeoModel("ModelA", "run1")
# unloadGeoModel("ModelA","run1")
# run_paths = runsInit("TestRun")
# createAstroParams(astroModels.bp_run1, ["p_temp","p_mag"])

# ORDER: Run initialization,

