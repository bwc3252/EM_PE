from __future__ import print_function
import numpy as np
from scipy.stats import multivariate_normal
import argparse

from em_pe.integrator_utils import monte_carlo_integrator

parser = argparse.ArgumentParser(description="Test monte carlo integrator")
parser.add_argument("--dim", type=int, default=1, help="Dimensions")
parser.add_argument("--cov", type=float, default=0.01, help="Cov")
parser.add_argument("--width", type=float, default=10.0, help="Width")
args = parser.parse_args()

def f(x):
    cov = args.cov * np.eye(args.dim)
    mean = np.zeros(args.dim)
    ret = np.rot90([multivariate_normal.logpdf(x, mean, cov)], -1)
    return ret


bounds = np.empty((args.dim, 2))
bounds[:,0] = -0.5 * args.width
bounds[:,1] = 0.5 * args.width
gmm_dict = {tuple(range(args.dim)):None}
ncomp = 1

integrator = monte_carlo_integrator.integrator(args.dim, bounds, gmm_dict, ncomp, use_lnL=True)
integrator.integrate(f, max_iter=40, progress=True, epoch=5)
print(integrator.integral, np.sqrt(integrator.var))
