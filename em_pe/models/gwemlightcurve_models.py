# -*- coding: utf-8 -*-
'''
Gwemlightcurves Models
----------------------
Models adapted from code found in `gwemlightcurves <https://github.com/mcoughlin/gwemlightcurves/tree/master/gwemlightcurves/KNModels/io>`_
'''
from __future__ import print_function
import numpy as np

from .model import model_base

class woko2017(model_base):
    '''
    Implementation of lightcurve model found `here <https://arxiv.org/abs/1705.07084>`_

    Parameters
    ----------
    weight : float
        Weight of the model
    '''

    def __init__(self, weight=1):
        name = 'woko2017'
        param_names = []
        bands = []
        model_base.__init__(self, name, param_names, bands, weight)

    def set_params(self, params, t_bounds):
        '''
        Method to set the parameters and run differential equation for lightcurve
        model.

        Parameters
        ----------
        params : dict
            Dictionary mapping parameter names to their values
        t_bounds : dict
            Dictionary mapping bands to [upper_bound, lower_bound] pairs for the
            time values observed in those bands
        '''
        pass


class me2017(model_base):
    '''
    Implementation of Metzger model

    Parameters
    ----------
    weight : float
        Weight of the model
    '''

    def __init__(self, weight=1):
        name = 'me2017'
        param_names = ['mej', 'vej']
        bands = []
        model_base.__init__(self, name, param_names, bands, weight)

    def set_params(self, params, t_bounds):
        '''
        Method to set the parameters and run differential equation for lightcurve
        model.

        Parameters
        ----------
        params : dict
            Dictionary mapping parameter names to their values
        t_bounds : dict
            Dictionary mapping bands to [upper_bound, lower_bound] pairs for the
            time values observed in those bands
        '''
