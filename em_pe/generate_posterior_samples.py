from __future__ import print_function
import numpy as np
import argparse

from models import model_dict, bounds_dict
from integrator_utils import monte_carlo_integrator

def parse_command_line_args():
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

def read_data(data_loc, files):
    data = {}
    for fname in files:
        band = fname.split('.')[0]
        data[band] = np.loadtxt(data_loc + fname)
    return data

def initialize_models(m):
    models = [model_dict[name](weight) for [name, weight] in m]
    ordered_params = []
    bands_used = []
    bounds = []
    for model in models:
        for param in model.param_names:
            if param not in ordered_params:
                ordered_params.append(param)
                bounds.append(bounds_dict[param])
        for band in model.bands:
            if band not in bands_used:
                bands_used.append(band)
    return models, ordered_params, bands_used, bounds

def evaluate_lnL(params, data, models, bands_used):
    temp_data = {}
    for band in bands_used:
        temp_data[band] = 0
    for model in models:
        model.set_params(params)
        for band in model.bands:
            t = data[band][0]
            temp_data[band] += model.evaluate(t, band)
    lnL = 0
    for band in bands_used:
        x = data[band][2]
        err = data[band][3]
        m = temp_data[band]
        diff = x - m
        lnL += np.sum(diff**2 / err**2)
    return -0.5 * lnL

def integrand(samples):
    n, _ = samples.shape
    ret = []
    for row in samples:
        params = dict(zip(ordered_params, row))
        ret.append(np.exp(evaluate_lnL(params, data, models, bands_used)))
    ret = np.array(ret).reshape((n, 1))
    return ret

def generate_samples(data, models, ordered_params, L_cutoff, bounds, min_iter, max_iter):
    # get a tuple of integers as dimensions to make the integrator happy
    param_ind = tuple(range(len(ordered_params)))
    dim = len(ordered_params) # number of dimensions
    k = 1 # number of gaussian components TODO make this something that can change
    gmm_dict = {param_ind:None}
    integrator = monte_carlo_integrator.integrator(dim, bounds, gmm_dict, k,
                    proc_count=None, L_cutoff=L_cutoff)
    integrator.integrate(integrand, min_iter=min_iter, max_iter=max_iter)
    samples = integrator.cumulative_values
    samples = np.append(samples, integrator.cumulative_p, axis=1)
    samples = np.append(samples, integrator.cumulative_p_s, axis=1)
    samples = np.append(samples, integrator.cumulative_samples, axis=1)
    return samples

if __name__ == '__main__':
    args = parse_command_line_args()
    m = args.m
    data_loc = args.dat
    files = args.f
    L_cutoff = args.cutoff
    min_iter = args.min
    max_iter = args.max
    if args.v:
        print('Loading data... ', end='')
    data = read_data(data_loc, files)
    if args.v:
        print('finished')
    models, ordered_params, bands_used, bounds = initialize_models(m)
    samples = generate_samples(data, models, ordered_params, L_cutoff, bounds, min_iter, max_iter)
    header = 'L p p_s ' + ' '.join(ordered_params)
    np.savetxt(args.out, samples, header=header)
