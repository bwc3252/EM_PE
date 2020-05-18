# -*- coding: utf-8 -*-
"""
Interpolated Model
------------------
Model interpolated from some surrogate model
"""
import numpy as np
import os

from .interpolated_model import interpolated

class me2017_lanthanide(interpolated):
    def __init__(self, weight=1):
        name = "me2017_lanthanide"
        param_names = ["mej", "vej"]
        bands = ["u", "g", "r", "i", "z", "y", "J", "H", "K"]
        interpolated.__init__(self, name, param_names, bands, weight)
