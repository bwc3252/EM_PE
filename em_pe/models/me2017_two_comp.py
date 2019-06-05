import numpy as np

from .model import model_base
from .me2017_lanthanide import me2017_lanthanide
from .me2017_non_lanthanide import me2017_non_lanthanide

class me2017_two_comp(model_base):
    def __init__(self, weight=1):
        name = "me2017_two_comp"
        param_names = ["mej_red", "mej_blue", "vej_red", "vej_blue", "frac"]
        bands = ["u", "g", "r", "i", "z", "y", "J", "H", "K"]
        model_base.__init__(self, name, param_names, bands, weight)
        self.red_model = me2017_lanthanide()
        self.blue_model = me2017_non_lanthanide()

    def set_params(self, params, t_bounds):
        self.params = params
        red_params = {"mej":params["mej_red"], "vej":params["vej_red"]}
        blue_params = {"mej":params["mej_blue"], "vej":params["vej_blue"]}
        self.red_model.set_params(red_params, t_bounds)
        self.blue_model.set_params(blue_params, t_bounds)

    def evaluate(self, tvec_days, band):
        red, _ = self.red_model.evaluate(tvec_days, band)
        blue, _ = self.blue_model.evaluate(tvec_days, band)
        frac = self.params["frac"]
        return (frac * red + (1.0 - frac) * blue), 0
