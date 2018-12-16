# -*- coding: utf-8 -*-
'''
Linear Model
------------
Simple model used for tests on actual data.
'''

from __future__ import print_function
import numpy as np
from .model import model_base

class linear_model(model_base):
    '''
    Simple model for tests on actual data

    Parameters
    ----------
    weight : float
        Weight of the model
    '''

    def __init__(self, weight=1):
        name = 'linear'
        param_names = ['m', 'y0']
        bands = ['R']
        model_base.__init__(self, name, param_names, bands, weight)

    def evaluate(self, t, band):
        '''
        Evaluate the model at specific time values using the current parameters.

        Parameters
        ----------
        t : np.ndarray
            Time values
        band : string
            Band to evaluate
        '''
        m = self.params['m']
        y0 = self.params['y0']
        # return lightcurve, model error
        return (m * t) + y0, np.zeros(len(t))
