# -*- coding: utf-8 -*-

import numpy as np
from scipy import interpolate
from gwemlightcurves.KNModels.io import DiUj2017

from .model import model_base

class diuj2017(model_base):
    def __init__(self, weight=1):
        name = "diuj2017"
        param_names = ["mej", "vej", "theta_ej", "phi_ej", "alpha"]
        bands = ["u", "g", "r", "i", "z", "y", "J", "H", "K"]
        model_base.__init__(self, name, param_names, bands, weight)
        wavelengths = np.array([354.3, 477.56, 612.95, 748.46, 865.78, 960.31,
            1235.0, 1662.0, 2159.0]) * 1.0e-7 # convert to cm (???)
        self.wavelength_dict = dict(zip(bands, wavelengths))
        self.index_dict = dict(zip(bands, range(len(bands))))
        self.mej = None
        self.vej = None
        self.theta_ej = None
        self.phi_ej = None
        self.alpha = None
    
    def set_params(self, params, t_bounds):
        self.mej = params["mej"]
        self.vej = params["vej"]
        self.theta_ej = params["theta_ej"]
        self.phi_ej = params["phi_ej"]
        self.alpha = params["alpha"]
    
    def evaluate(self, tvec_days, band):
        tmin = 0.1
        tmax = tvec_days[-1]
        dt = 0.05

        ### hard-code some parameters
        vmin = 0.02
        kappa = 1.0
        eps = 1.58e10
        eth = 0.5
        flgbct = 1.0
        beta = 3.0
        slope_r = -1.2
        theta_r = 0.0
        
        tdays, L, lightcurves = DiUj2017.calc_lc(tmin, tmax, dt, self.mej,
                self.vej, vmin, self.theta_ej, self.phi_ej,
                kappa, eps, self.alpha, eth, flgbct)
        
        lc = lightcurves[self.index_dict[band]]
        mask = np.isfinite(lc) # filter out NaNs, since interp1d is undefined for them
        f = interpolate.interp1d(tdays[mask], lc[mask], fill_value="extrapolate")
        return f(tvec_days), 0
