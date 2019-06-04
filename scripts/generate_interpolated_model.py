# -*- coding: utf-8 -*-
"""
Generate Interpolated Models
----------------------------
"""
from __future__ import print_function
import argparse
import os
import numpy as np
from sklearn.gaussian_process import GaussianProcessRegressor
from joblib import dump
import pyDOE
from scipy import interpolate

### surrogate model imports
from gwemlightcurves.KNModels.io import Me2017

### parse command line arguments
parser = argparse.ArgumentParser(description="Generate a model file for an interpolated model")
parser.add_argument("--m", default="", help="Name of model")
parser.add_argument("--tmin", type=float, default=0.0, help="Start time for interpolated model")
parser.add_argument("--tmax", type=float, default=10.0, help="Stop time for interpolated model")
parser.add_argument("--n", type=int, default=1000, help="Number of samples to use")
args = parser.parse_args()

class surrogate_model:
    def __init__(self, args):
        self.args = args
        self.name = args.m

        self.lc_func = None
        self.param_names = None
        self.lims = None
        self.bands = None
        self.kappa_r = None

        ### check if available model is being requested
        if self.name == "me2017_non_lanthanide":
            self.lc_func = self.me2017
            self.param_names = ["mej", "vej"]#, "kappa_r_nl"]
            self.lims = np.array([[1.0e-7, 0.1],            # mej
                                  [0.01, 0.4],              # vej
                                 # [0.7, 5.0],               # kappa_r
                                  [args.tmin, args.tmax]])  # t
            self.kappa_r = 1.0
            self.bands = ["u", "g", "r", "i", "z", "y", "J", "H", "K"]
        elif self.name == "me2017_lanthanide":
            self.lc_func = self.me2017
            self.param_names = ["mej", "vej"]#, "kappa_r_l"]
            self.lims = np.array([[1.0e-7, 0.1],            # mej
                                  [0.01, 0.4],              # vej
                                 # [5.0, 15.0],              # kappa_r
                                  [args.tmin, args.tmax]])  # t
            self.kappa_r = 10.0
            self.bands = ["u", "g", "r", "i", "z", "y", "J", "H", "K"]
        else:
            print("Error: model '" + self.name + "' does not exist")
            exit()

    def me2017(self, band, mej, vej, t):
        dt = 0.05
        ret = np.empty(args.n)
        index_dict = dict(zip(self.bands, range(len(self.bands))))
        beta = 3.0
        kappa_r = self.kappa_r
        for i in range(args.n):
            _mej, _vej, _t = mej[i], vej[i], t[i]
            tdays, L, lightcurves, Tobs = Me2017.calc_lc(args.tmin, args.tmax, dt, _mej, _vej, beta, kappa_r)
            lc = lightcurves[index_dict[band]]
            mask = np.isfinite(lc) # filter out NaNs, since interp1d is undefined for them
            f = interpolate.interp1d(tdays[mask], lc[mask], fill_value="extrapolate")
            ret[i] = f(_t)
        return ret


### initialize model
model = surrogate_model(args)

ndim, _ = model.lims.shape

### generate Latin hypercube samples and transform to correct intervals
x = pyDOE.lhs(ndim, samples=args.n) - 0.5
interval_lengths = model.lims[:,1] - model.lims[:,0]
x *= interval_lengths
centers = np.mean(model.lims, axis=1).flatten()
x += centers

for band in model.bands:
    print("Sampling and interpolating", band, "band...")
    ### evaluate function at randomly sampled points
    y = model.lc_func(band, *x.T)
    
    ### initialize model
    gp = GaussianProcessRegressor()

    ### fit model
    gp.fit(x, y)

    ### write model to disk
    fname = os.environ["EM_PE_INSTALL_DIR"] + "/Data/" + args.m + "_" + band + "_model.joblib"
    dump(gp, fname)
