# -*- coding: utf-8 -*-
'''
Oriented model
--------------
'''

import numpy as np

from . import *
 

class oriented(model_base):
    '''
    Model class to wrap another model and provide orientation dependance
    functionality
    '''
    def __init__(self, model_name, profile='flat', weight=1):
        self.model = model_dict[model_name](weight)
        self.profile = profile
        self.weight = weight
        name = 'oriented'
        param_names = self.model.param_names + ['angle']
        bands = self.model.bands
        model_base.__init__(self, name, param_names, bands, weight)
        self.params = None
        self.max_angle = 10 # degrees

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
        self.model.set_params(params, t_bounds)
    
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
        m, err = self.model.evaluate(tvec_days, band)
        if self.profile == 'flat':
            if self.params['angle'] > self.max_angle:
                m = np.zeros(m.size)
        elif self.profile == 'gaussian':
            m *= np.exp(-1 * self.params['angle']**2)
        return m, err

