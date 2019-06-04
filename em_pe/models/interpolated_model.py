# -*- coding: utf-8 -*-
"""
Interpolated Model
------------------
Model interpolated from some surrogate model
"""
from __future__ import print_function
import numpy as np
import os
from joblib import load

from .model import model_base

class interpolated(model_base):
    def __init__(self, name, param_names, bands, weight=1):
        model_base.__init__(self, name, param_names, bands, weight)
        self.gp_dict = None

    def evaluate(self, tvec_days, band):
        gp = self.gp_dict[band]
        x = np.empty((len(tvec_days), len(self.param_names) + 1))
        for i in range(len(self.param_names)):
            x[:,i] = self.params[self.param_names[i]]
        x[:,i + 1] = tvec_days
        return gp.predict(x, return_std=True)
