# -*- coding: utf-8 -*-
"""
Lightcurve Utils
----------------
Functions for converting BNS parameters to ejecta parameters

Mostly taken from `here <https://github.com/mcoughlin/gwemlightcurves/blob/master/gwemlightcurves/EjectaFits/Di2018b.py>`_
"""

import numpy as np
import scipy

def calc_mej(m1, lambda1, m2, lambda2):
    """
    Funtion to calculate ejecta mass based on BNS parameters
    """
    c1, c2 = calc_compactness(lambda1, lambda2)

    a= -0.0719
    b= 0.2116
    d= -2.42
    n= -2.905

    log10_mej = a*(m1*(1-2*c1)/c1 + m2*(1-2*c2)/c2) + b*(m1*(m2/m1)**n + m2*(m1/m2)**n)+d
    meje_fit = 10**log10_mej
    return meje_fit	

def calc_vej(m1, lambda1, m2, lambda2):
    """
    Funtion to calculate ejecta velocity based on BNS parameters
    """
    c1, c2 = calc_compactness(lambda1, lambda2)

    a=-0.3090
    b=0.657
    c=-1.879
    return a*(m1/m2)*(1+c*c1) + a*(m2/m1)*(1+c*c2)+b

def calc_compactness(lambda1, lambda2):
    """
    Funtion to calculate compactness based on lambda1 and lambda2
    """
    ### calculate compacness from parameters: see https://arxiv.org/pdf/1812.04803.pdf, Appendix B II
    c1 = 0.371 - 0.0391 * np.log(lambda1) + 0.001056 * np.log(lambda1)**2
    c2 = 0.371 - 0.0391 * np.log(lambda2) + 0.001056 * np.log(lambda2)**2

    return c1, c2
