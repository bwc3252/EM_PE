from __future__ import print_function
import numpy as np

from .model import model_base

class one_band_test_model(model_base):

    def __init__(self):
        name = 'one_band_test_model'
        param_names = ['a', 'b']
        bands = ['test_bandA']
        model_base.__init__(self, name, param_names, bands)

    def evaluate(self, t):
        a = self.params['a']
        b = self.params['b']
        vals = (1 / np.cosh(a * t)) + b
        return {'test_bandA':vals}


class two_band_test_model(model_base):

    def __init__(self):
        name = 'two_band_test_model'
        param_names = ['a', 'b']
        bands = ['test_bandA', 'test_bandB']
        model_base.__init__(self, name, param_names, bands)

    def evaluate(self, t):
        a = self.params['a']
        b = self.params['b']
        A_vals = (1 / np.cosh(a * t)) + b
        B_vals = (1 / np.cosh(a * t)) - b
        return {'test_bandA':A_vals, 'test_bandB':B_vals}
