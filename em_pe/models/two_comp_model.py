# -*- coding: utf-8 -*-
'''
Two-Component model
-------------------
'''

import numpy as np

from .model import model_base
from .woko2017 import woko2017

class two_comp(model_base):

    def __init__(self, weight=1):
        name = 'two_comp'
        param_names = ['mej1', 'vej1', 'mej2', 'vej2', 'frac', 'dist']
        bands = ['g', 'r', 'i', 'z', 'y', 'J', 'H', 'K']
        model_base.__init__(self, name, param_names, bands, weight)
        self.model1 = woko2017(1.0)
        self.model2 = woko2017(1.0)

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
        model1_params = {'mej':params['mej1'], 'vej':params['vej1']}
        model2_params = {'mej':params['mej2'], 'vej':params['vej2']}
        self.model1.set_params(model1_params, t_bounds)
        self.model2.set_params(model2_params, t_bounds)

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
        m1, m1err = self.model1.evaluate(tvec_days, band)
        m2, m2err = self.model2.evaluate(tvec_days, band)
        w1 = self.params['frac']
        w2 = 1 - w1
        return ((w1 * m1) + (w2 * m2)), 0
