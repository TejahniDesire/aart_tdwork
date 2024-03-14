import subprocess
import EZPaths
import os
from aart_func import *
from params import *
import astroModels
import fileloading

all_models, string = fileloading.createAstroParams(astroModels.bp_run1, ["p_temp","p_mag"])
print(string)