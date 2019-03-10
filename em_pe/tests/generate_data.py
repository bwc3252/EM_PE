from __future__ import print_function
import numpy as np
import argparse
import sys

from em_pe.models import model_dict

parser = argparse.ArgumentParser(description='Generate synthetic data for PE tests')
parser.add_argument('--m', help='Name of model to use')
parser.add_argument('--out', help='Directory to output data files to')
parser.add_argument('--p', action='append', nargs=2, help='Parameter value. Must set each parameter value for the model being used, as a pair (e.g. "--p dist 40")')
parser.add_argument('--tmin', type=float, help='Minimum time (in days)')
parser.add_argument('--tmax', type=float, help='Maximum time (in days)')
parser.add_argument('--n', type=int, help='Number of data points per band')
parser.add_argument('--err', type=float, help='Error std. dev.')
args = parser.parse_args()

### initialize the model
model = model_dict[args.m]()

### create a dictionary mapping parameter names to values
params = dict(args.p)
for p in params:
    params[p] = float(params[p])

t_bounds = [args.tmin, args.tmax]

### set model parameters
model.set_params(params, t_bounds)

### generate times
tdays = np.linspace(args.tmin, args.tmax, args.n)

### generate and save the data
for band in model.bands:
    filename = args.out + band + '.txt'
    m, _ = model.evaluate(tdays, band)
    data = np.empty((4, args.n))
    data[0] = tdays
    data[2] = m + np.random.uniform(-1 * args.err, args.err, args.n) # generate errors
    data[3] = np.ones(args.n) * args.err
    np.savetxt(filename, data.T)

### save true values

truths = []
for param in model.param_names:
    truths.append(params[param])
filename = args.out + 'test_truths.txt'
np.savetxt(filename, truths)
