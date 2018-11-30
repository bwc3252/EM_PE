from __future__ import print_function
import numpy as np

from .model import model_base

class one_band_test_model(model_base):
    '''
    Simple one-band model for testing
    '''

    def __init__(self, weight=1):
        name = 'one_band_test_model'
        param_names = ['a', 'b']
        bands = ['test_bandA']
        model_base.__init__(self, name, param_names, bands, weight)

    def evaluate(self, t, band):
        a = self.params['a']
        b = self.params['b']
        ret = (1 / np.cosh(a * t)) + b
        return ret #* self.weight


class two_band_test_model(model_base):
    '''
    Simple two-band model for testing
    '''

    def __init__(self, weight=1):
        name = 'two_band_test_model'
        param_names = ['a', 'b']
        bands = ['test_bandA', 'test_bandB']
        model_base.__init__(self, name, param_names, bands, weight)

    def evaluate(self, t, band):
        a = self.params['a']
        b = self.params['b']
        if band == 'test_bandA':
            ret = (1 / np.cosh(a * t)) + b
        else:
            ret = (1 / np.cosh(a * t)) - b
        return ret * self.weight
