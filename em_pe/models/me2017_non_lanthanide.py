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

from .interpolated_model import interpolated

class me2017_non_lanthanide(interpolated):
    def __init__(self, weight=1):
        name = "me2017_non_lanthanide"
        param_names = ["mej", "vej"]#, "kappa_r_nl"]
        bands = ["u", "g", "r", "i", "z", "y", "J", "H", "K"]
        interpolated.__init__(self, name, param_names, bands, weight)
        self.gp_dict = {}
        self.fname_base = os.environ["EM_PE_INSTALL_DIR"] + "/Data/me2017_non_lanthanide_" 
        for band in bands:
            fname = self.fname_base + band + "_model.joblib"
            self.gp_dict[band] = load(fname)
