import numpy as np

from .model import model_base

class test(model_base):
    def __init__(self, weight=1):
        name = "test"
        bands = ["test_band"]
        param_names = ["a", "b", "c"]
        model_base.__init__(self, name, param_names, bands, weight)
    
    def evaluate(self, tvec_days, band):
        a = self.params["a"]
        b = self.params["b"]
        c = self.params["c"]
        return a * np.log(b * tvec_days + 1) + c, 0
