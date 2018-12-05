# -*- coding: utf-8 -*-
'''
Model
-----
Base class for lightcurve models
'''
from __future__ import print_function

class model_base:
    '''
    Base "template" class for lightcurve models

    Parameters
    ----------
    name : string
        Name of the model
    param_names : list
        Names of parameters
    bands : list
        Names of data bands
    weight : float
        Weight of model used in likelihood calculations
    '''
    def __init__(self, name, param_names, bands, weight=1):
        self.name = name
        self.param_names = param_names
        self.bands = bands
        self.weight = weight
        self.params = None
        self.t_bounds = None

    def set_params(self, params, t_bounds):
        '''
        Method to set the parameters. This would be overridden for more complex
        models, because this is where the differential equations would be solved.

        Parameters
        ----------
        params : dict
            Dictionary mapping parameter names to their values
        t_bounds : dict
            Dictionary mapping bands to [upper_bound, lower_bound] pairs for the
            time values observed in those bands
        '''
        self.params = params
        self.t_bounds = t_bounds

    ### Functions that should be implemented by child classes.
    ### NOTE: These could be seen as (and maybe should be) abstract methods

    def evaluate(self, t, band):
        '''
        Method to evaluate model at specific time values using the current parameters.
        Must be implemented by all child classes.

        Parameters
        ----------
        t : np.ndarray
            Time values
        band : string
            Band to evaluate
        '''
        pass
