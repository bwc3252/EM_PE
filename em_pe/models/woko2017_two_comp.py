# -*- coding: utf-8 -*-
'''
Two-Component model
-------------------
'''

import numpy as np

from .model import model_base
from .woko2017 import woko2017

class woko2017_two_comp(model_base):

    def __init__(self, weight=1):
        name = 'two_comp'
        param_names = ['mej_red', 'vej_red', 'mej_blue', 'vej_blue', 'frac']
        bands = ['g', 'r', 'i', 'z', 'y', 'J', 'H', 'K']
        model_base.__init__(self, name, param_names, bands, weight)
        self.blue_model = woko2017(kappa_r=1.0)
        self.red_model = woko2017(kappa_r=10.0)

    def set_params(self, params, t_bounds):
        '''
        Method to set the parameters for lightcurve
        model.

        Parameters
        ----------
        params : dict
            Dictionary mapping parameter names to their values
        t_bounds : list
            [upper bound, lower bound] pair for time values
        '''
        self.params = params
        blue_params = {'mej':params['mej_blue'], 'vej':params['vej_blue']}
        red_params = {'mej':params['mej_red'], 'vej':params['vej_red']}
        self.blue_model.set_params(blue_params, t_bounds)
        self.red_model.set_params(red_params, t_bounds)

    def evaluate(self, tvec_days, band):
        '''
        Evaluate model at specific time values using the current parameters.

        Parameters
        ----------
        tvec_days : np.ndarray
            Time values
        band : string
            Band to evaluate
        '''
        blue, _ = self.blue_model.evaluate(tvec_days, band)
        red, _ = self.red_model.evaluate(tvec_days, band)
        w1 = self.params['frac']
        w2 = 1 - w1
        return ((w1 * red) + (w2 * blue)), 0
