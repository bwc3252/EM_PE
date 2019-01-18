# -*- coding: utf-8 -*-
'''
Generate Posterior Samples
--------------------------
Generate posterior parameter samples according to the given command line arguments

Example::

    $ python generate_posterior_samples.py --dat Data/ --m test_modelA --cutoff 10e-8 --f test_bandA.txt --out posterior_samples.txt

To see full command line parameter documentation::

    $ python generate_posterior_samples.py -h
    usage: generate_posterior_samples.py [-h] [--dat DAT] [--m M M] [-v]
                                     [--cutoff CUTOFF] [--f F] [--min MIN]
                                     [--max MAX] [--out OUT]

    Generate posterior parameter samples from lightcurve data

    optional arguments:
      -h, --help       show this help message and exit
      --dat DAT        Path to data directory
      --m M M          Name of a model to use
      -v               Verbose mode
      --cutoff CUTOFF  Likelihood cutoff for storing posterior samples
      --f F            Name of a data file
      --min MIN        Minimum number of integrator iterations
      --max MAX        Maximum number of integrator iterations
      --out OUT        Location to store posterior samples
'''
from __future__ import print_function
import numpy as np
import argparse
import sys

from models import model_dict, bounds_dict
from integrator_utils import monte_carlo_integrator

def _parse_command_line_args():
    '''
    Parses and returns the command line arguments.
    '''
    parser = argparse.ArgumentParser(description='Generate posterior parameter samples from lightcurve data')
    parser.add_argument('--dat', help='Path to data directory')
    parser.add_argument('--m', action='append', nargs=2, help='Name of a model to use')
    parser.add_argument('-v', action='store_true', help='Verbose mode')
    parser.add_argument('--cutoff', default=0, type=float, help='Likelihood cutoff for storing posterior samples')
    parser.add_argument('--f', action='append', help='Name of a data file')
    parser.add_argument('--min', default=20, type=int, help='Minimum number of integrator iterations')
    parser.add_argument('--max', default=20, type=int, help='Maximum number of integrator iterations')
    parser.add_argument('--out', help='Location to store posterior samples')
    return parser.parse_args()

def _read_data(data_loc, files):
    if args.v:
        print('Loading data... ', end='')
    data = {}
    bands_used = []
    for fname in files:
        band = fname.split('.txt')[0]
        data[band] = np.loadtxt(data_loc + fname)
        bands_used.append(band)
    if args.v:
        print('finished')
    return data, bands_used

def _initialize_models(m, bands_used):
    if args.v:
        print('Initializing models... ', end='')
    ### initialize model objects
    models = [model_dict[name](float(weight)) for [name, weight] in m]
    ordered_params = [] # keep track of all parameters used
    bounds = [] # bounds for each parameter
    for model in models:
        for param in model.param_names:
            if param not in ordered_params:
                ordered_params.append(param)
                bounds.append(bounds_dict[param])
    t_bounds = [np.inf, -1 * np.inf] # tmin and tmax for each band
    for band in bands_used:
        t = data[band][0]
        t_bounds[0] = min(min(t), t_bounds[0])
        t_bounds[1] = max(max(t), t_bounds[1])
    if args.v:
        print('finished')
    return models, ordered_params, bounds, t_bounds

def _evaluate_lnL(params, data, models, bands_used, t_bounds):
    temp_data = {} # used to hold model data and squared model error
    for band in bands_used:
        temp_data[band] = [0, 0]
    ### evaluate each model for the params, in the required bands
    for model in models:
        model.set_params(params, t_bounds)
        for band in model.bands:
            if band in bands_used:
                t = data[band][0]
                m, m_err = model.evaluate(t, band)
                temp_data[band][0] += m * model.weight
                temp_data[band][1] += m_err**2
    lnL = 0
    ### calculate lnL
    for band in bands_used:
        x = data[band][2]
        err = data[band][3]
        m = temp_data[band][0]
        m_err_squared = temp_data[band][1]
        diff = x - m
        lnL += np.sum(diff**2 / (err**2 + m_err_squared))
    return -0.5 * lnL

def _integrand(samples):
    n, _ = samples.shape
    ret = []
    for row in samples:
        params = dict(zip(ordered_params, row)) # map each parameter's name to its value
        ret.append(_evaluate_lnL(params, data, models, bands_used, t_bounds))
    ret = np.array(ret).reshape((n, 1))
    return ret

def generate_samples(data, models, ordered_params, L_cutoff, bounds, min_iter, max_iter):
    '''
    Generate posterior samples

    Parameters
    ----------
    data : dict
        Dictionary mapping bands to data
    models : dict
        Dictionary mapping model names to objects
    ordered_params : list
        List of parameter names
    L_cutoff : float
        Likelihood cutoff for storing samples
    bounds : np.ndarray
        Limits of integration
    min_iter : int
        Minimum number of integrator iterations
    max_iter : int
        Maximum number of integrator iterations
    '''
    if args.v:
        print('Generating posterior samples')
        sys.stdout.flush()
    ### get a tuple of integers as dimensions to make the integrator happy
    param_ind = tuple(range(len(ordered_params)))
    dim = len(ordered_params) # number of dimensions
    k = 1 # number of gaussian components TODO make this something that can change
    gmm_dict = {param_ind:None}
    ### initialize and run the integrator
    integrator = monte_carlo_integrator.integrator(dim, bounds, gmm_dict, k,
                    proc_count=None, L_cutoff=L_cutoff, use_lnL=True)
    integrator.integrate(_integrand, min_iter=min_iter, max_iter=max_iter)
    ### make the array of samples
    samples = integrator.cumulative_values
    samples = np.append(samples, integrator.cumulative_p, axis=1)
    samples = np.append(samples, integrator.cumulative_p_s, axis=1)
    samples = np.append(samples, integrator.cumulative_samples, axis=1)
    if args.v:
        print('Integral result:', integrator.integral)
    return samples

if __name__ == '__main__':
    args = _parse_command_line_args()
    m = args.m
    data_loc = args.dat
    files = args.f
    L_cutoff = args.cutoff
    min_iter = args.min
    max_iter = args.max
    data, bands_used = _read_data(data_loc, files)
    models, ordered_params, bounds, t_bounds = _initialize_models(m, bands_used)
    samples = generate_samples(data, models, ordered_params, L_cutoff, bounds, min_iter, max_iter)
    header = 'lnL p p_s ' + ' '.join(ordered_params)
    np.savetxt(args.out, samples, header=header)
