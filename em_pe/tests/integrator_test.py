from __future__ import print_function
import numpy as np
from scipy.stats import multivariate_normal

from em_pe.integrator_utils import monte_carlo_integrator

def f(x):
    
    cov = np.array([[0.01, 0, 0],
                    [0, 0.01, 0],
                    [0, 0, 0.01]])
    mean = np.array([0, 0, 0])
    ret = np.rot90([multivariate_normal.logpdf(x, mean, cov)], -1)
    return ret #np.ones((len(x), 1))


d = 3 # three dimensions
bounds = np.array([[-1.0, 1.0],
                   [-1.0, 1.0],
                   [-1.0, 1.0]])
gmm_dict = {(0, 1, 2):None}
ncomp = 1
integrator = monte_carlo_integrator.integrator(d, bounds, gmm_dict, ncomp, use_lnL=True)
integrator.integrate(f, max_iter=40, progress=True)
print(integrator.integral, np.sqrt(integrator.var))
