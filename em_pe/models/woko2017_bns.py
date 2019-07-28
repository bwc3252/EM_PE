# -*- coding: utf-8 -*-
'''
Woko2017 Model
----------------------
Models adapted from code found in `gwemlightcurves <https://github.com/mcoughlin/gwemlightcurves/tree/master/gwemlightcurves/KNModels/io>`_
'''
from __future__ import print_function
import numpy as np
from scipy.interpolate import splrep, splev
import os

from .model import model_base
from .woko2017 import woko2017

from em_pe.utils import calc_mej, calc_vej, calc_compactness
import EOSManager

class woko2017_bns(model_base):
    '''
    Implementation of lightcurve model found `here <https://arxiv.org/abs/1705.07084>`_

    Parameters
    ----------
    weight : float
        Weight of the model
    '''

    def __init__(self, weight=1, kappa_r=10.0):
        name = "woko2017_bns"
        param_names = ["m1", "m2"]
        bands = ["g", "r", "i", "z", "y", "J", "H", "K"]
        model_base.__init__(self, name, param_names, bands, weight)
        self.base_model = woko2017()
        self.eos = EOSManager.EOSLALSimulation(name="AP4")
    
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
        m1, m2 = params["m1"], params["m2"]
        lambda1 = self.eos.lambda_from_m(params["m1"])
        lambda2 = self.eos.lambda_from_m(params["m2"])
        mej = calc_mej(m1, lambda1, m2, lambda2)
        vej = calc_vej(m1, lambda1, m2, lambda2)
        self.base_model.set_params({"mej":mej, "vej":vej}, t_bounds)

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
        return self.base_model.evaluate(tvec_days, band)
