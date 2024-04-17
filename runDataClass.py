import astroModels


class runData():

    def __init__(self, name:str ,brightparams_grid: dict,variable_brightparams_names:list,
                 variable_geo_values,variable_geo_names):
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

    def setAction(self,action):
        """

        Args:
            action: {
                "var": str,
                "start": Float,
                "stop": Float,
                "step": Float,
                "images": Bolean

        """
        self.action=action

    def setisNormalized(self,isnormalized):
        self.isNormalized = isnormalized

    def setBlurrPolicy(self,blurr_policy):
        self.blurr_policy = blurr_policy

    def getRunName(self):
        return self.run_name

    def getBrightparams(self):
        return self.brightparams_values_grid

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
               [(["a"], [str(.3)]),(["a"], [str(.9)])],
               ["ModelA", "ModelB"])
# testRun1______________________________________________________________________

testRun1 = runData("testRun1",
               astroModels.bp_testRun1,
               ["p_temp", "p_mag"],
               [(["a"], [str(.3)]),(["a"], [str(.9)])],
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
                   [(["a"], [str(.3)]),(["a"], [str(.9)])],
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
                   [(["a"], [str(.3)]),(["a"], [str(.9)])],
                   ["ModelA", "ModelB"])

# soloRun2______________________________________________________________________
soloRun2 = runData("soloRun2",
                   astroModels.bp_soloRun2,
                   [],
                   [(["a"], [str(.9)])],
                   ["ModelB"])
testRun2.setAction(
    {
            "var": "nu0",
            "start": 670e9,
            "stop": 700e9,
            "step": 20e9,
            "images": True
    }
)

