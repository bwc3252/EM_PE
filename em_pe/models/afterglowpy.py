# -*- coding: utf-8 -*-
"""
Interpolated Model
------------------
Model interpolated from some surrogate model
"""
from __future__ import print_function
import numpy as np
import os

from .interpolated_model import interpolated

class afterglowpy(interpolated):
    def __init__(self, weight=1):
        name = "afterglowpy"
        param_names = ["log_E0", "thetaV"]
        bands = ["u", "g", "r", "i", "z", "y", "J", "H", "K"]
        interpolated.__init__(self, name, param_names, bands, weight)
