# -*- coding: utf-8 -*-
"""
Validate Interpolated Models
----------------------------
"""
from __future__ import print_function
import argparse
import os
import numpy as np
from scipy import interpolate

from em_pe.models import model_dict, bounds_dict

parser = argparse.ArgumentParser(description="Check the accuracy of a previously interpolated model")
parser.add_argument("--m", default="", help="Name of model")
parser.add_argument("--tmin", type=float, default=0.0, help="Start time for interpolated model")
parser.add_argument("--tmax", type=float, default=10.0, help="Stop time for interpolated model")
parser.add_argument("--n", type=int, default=1000, help="Number of samples to use")
parser.add_argument("--n-times", type=int, default=50, help="Number of time points to use")
args = parser.parse_args()

if args.m in ["me2017_lanthanide", "me2017_non_lanthanide"]:
    from gwemlightcurves.KNModels.io import Me2017
elif args.m in ["afterglowpy"]:
    import grbpy as grb

class surrogate_model:
    def __init__(self, args):
        self.args = args
        self.name = args.m

        self.lc_func = None
        self.param_names = None
        self.lims = None
        self.bands = None
        self.fixed_params ={}

        ### check if available model is being requested
        if self.name == "me2017_non_lanthanide":
            self.lc_func = self.me2017
            self.param_names = ["mej", "vej"]
            self.lims = np.array([[1.0e-5, 0.5],            # mej
                                  [0.001, 0.9]])            # vej
            self.fixed_params["kappa_r"] = 1.0
            self.bands = ["u", "g", "r", "i", "z", "y", "J", "H", "K"]
        elif self.name == "me2017_lanthanide":
            self.lc_func = self.me2017
            self.param_names = ["mej", "vej"]
            self.lims = np.array([[1.0e-5, 0.5],            # mej
                                  [0.001, 0.9]])            # vej
            self.fixed_params["kappa_r"] = 10.0
            self.bands = ["u", "g", "r", "i", "z", "y", "J", "H", "K"]
        elif self.name == "afterglowpy":
            self.lc_func = self.afterglowpy
            self.param_names = ["log_E0", "thetaV"]
            self.lims = np.array([[70.0, 125.0],            # log_E0
                                  [0.0, np.pi / 2]])        # thetaV
            self.bands = ["u", "g", "r", "i", "z", "y", "J", "H", "K"]
        else:
            print("Error: model '" + self.name + "' does not exist")
            exit()

    def me2017(self, band, t, mej, vej):
        dt = 0.05
        index_dict = dict(zip(self.bands, range(len(self.bands))))
        beta = 3.0
        kappa_r = self.fixed_params["kappa_r"]
        tdays, L, lightcurves, Tobs = Me2017.calc_lc(args.tmin, args.tmax, dt, 
                mej, vej, beta, kappa_r)
        lc = lightcurves[index_dict[band]]
        mask = np.isfinite(lc) # filter out NaNs, since interp1d is undefined for them
        f = interpolate.interp1d(tdays[mask], lc[mask], fill_value="extrapolate")
        return f(t)

    def afterglowpy(self, band, t, log_E0, thetaV):
        ### constants
        day = 86400.0
        pc = 3.086e18
        ### hard-coded parameters (see afterglowpy for more information)
        jetType = 0     # for "Gaussian" shape
        specType = 0
        thC = 0.08
        thW = 0.35
        b = 6
        L0 = 0.0
        q = 0.0
        ts = 0.0
        n0 = 1.0e-3
        p = 2.15
        epse = 1.0e-1
        epsB = 1.0e-3
        ksiN = 1.0
        dL = 10 * pc    # fiducial distance, since we'll add the distance modulus later
        ### get lightcurve
        t_sec = t * day
        E0 = np.exp(log_E0)
        Y = np.array([thetaV, E0, thC, thW, b, L0, q, ts, n0, p, epse, epsB, ksiN, dL])
        wavelengths = np.array([354.3, 477.56, 612.95, 748.46, 865.78, 960.31, 1235.0, 1662.0, 2159.0]) * 1.0e-7 # convert to cm (???)
        wavelength_dict = dict(zip(self.bands, wavelengths))
        nu = np.ones(t.size) * 3.0e10 / wavelength_dict[band]
        Fnu = grb.fluxDensity(t_sec, nu, jetType, specType, *Y, spread=True, latRes=5)
        mAB = -2.5 * np.log10(Fnu) - 48.60
        return mAB

m_interp = model_dict[args.m]()
m_surr = surrogate_model(args)

bounds = np.array([bounds_dict[p] for p in m_interp.param_names])
params = np.random.uniform(bounds[:,0], bounds[:,1], (args.n, len(m_interp.param_names)))

s = np.zeros(args.n)
for i in range(args.n):
    p = params[i]
    t = np.sort(np.random.uniform(args.tmin, args.tmax, args.n_times))
    p_dict = dict(zip(m_surr.param_names, p))
    m_interp.set_params(p_dict, [args.tmin, args.tmax])
    b = np.random.choice(m_interp.bands)
    diff = m_surr.lc_func(b, t, *p) - m_interp.evaluate(t, b)[0]
    s[i] += np.sum(diff**2)

print("median summed squared difference:", np.median(s))
print("max summed squared difference:", np.max(s))
