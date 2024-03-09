import subprocess
import EZPaths
import os
from aart_func import *
from params import *


def intensityNameWrite(brightparams,funckeys):
    filename = path + ('Intensity_a_{}_i_{}_nu_{}_mass_{}_scaleh_{}_thetab_{}_beta_{}_rie_{}_rb_{}_nth0_{}_te0_{}_'
                       'b0_{}_pdens_{}_ptemp_{}_pmag_{}_nscale_{}_emkey_{}_bkey_{}_nkey_{}_tnkey_{}_bnkey_{}.h5').format(
        spin_case,
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


def intensityNameRead(brightparams,funckeys):

    filename = path + ('Intensity_a_{}_i_{}_nu_{}_mass_{}_scaleh_{}_thetab_{}_beta_{}_rie_{}_rb_{}_nth0_{}_te0_{}_'
                       'b0_{}_pdens_{}_ptemp_{}_pmag_{}_nscale_{}_emkey_{}_bkey_{}_nkey_{}_tnkey_{}_bnkey_{}').format(
        spin_case,
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


def loadGeoModel(model_name:str, run:str):
    paramFile = "./params.py"
    model_is_loaded = os.path.exists(paramFile)
    if model_is_loaded:
        print("MODEL ALREADY LOADED, REMOVING PREVIOUS MODEL")
        subprocess.run(["rm " + paramFile], shell=True)

    cmd = "cp " + EZPaths.modeldir + run + "/geoParams/" + model_name + ".py" + " ./params.py"

    subprocess.run([cmd], shell=True)


def unloadGeoModel(model_name:str, run):

    cmd = "mv " + "./params.py" + " " + EZPaths.modeldir + run + "/geoParams/" + model_name + ".py"
    print(cmd)
    subprocess.run([cmd], shell=True)




#
# loadGeoModel("ModelA", "run1")
# unloadGeoModel("ModelA","run1")

