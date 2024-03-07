import subprocess
import EZPaths
import os


def loadModel(model_name:str, run:str):
    paramFile = "./params.py"
    model_is_loaded = os.path.exists(paramFile)
    if model_is_loaded:
        print("MODEL ALREADY LOADED, REMOVING PREVIOUS MODEL")
        subprocess.run(["rm " + paramFile], shell=True)

    cmd = "cp " + EZPaths.modeldir + run + "/" + model_name + ".py" + " ./params.py"

    subprocess.run([cmd], shell=True)


def unloadModel(model_name:str, run):

    cmd = "mv " + "./params.py" + " " + EZPaths.modeldir + run + "/" + model_name + ".py"
    print(cmd)
    subprocess.run([cmd], shell=True)

#
# loadModel("ModelA", "run1")
unloadModel("ModelA","run1")

