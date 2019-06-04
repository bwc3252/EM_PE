# -*- coding: utf-8 -*-
'''
Generate test data from a model
-------------------------------
Code to generate synthetic test data from one of the models.
'''

from __future__ import print_function
import numpy as np
import argparse
import sys
from astropy.time import Time
import json

import lal
import lalsimulation as lalsim

from em_pe.models import model_dict
from em_pe.lightcurve_utils import calc_mej, calc_vej

parser = argparse.ArgumentParser(description='Generate synthetic data for PE tests')
parser.add_argument('--m', help='Name of model to use')
parser.add_argument('--out', help='Directory to output data files to')
parser.add_argument('--p', action='append', nargs=2, help='Parameter value. Must set each parameter value for the model being used, as a pair (e.g. "--p dist 40")')
parser.add_argument('--tmin', type=float, help='Minimum time (in days)')
parser.add_argument('--tmax', type=float, help='Maximum time (in days)')
parser.add_argument('--n', type=int, help='Number of data points per band')
parser.add_argument('--err', type=float, help='Error std. dev.')
parser.add_argument('--orientation', nargs=2, help='Orientation profile and angle')
parser.add_argument('--t0', type=float, default=0, help='Start time for event')
parser.add_argument('--time_format', default='gps', help='Time format for t0 (gps or mjd)')
parser.add_argument('--json', action='store_true', help='Export data in JSON format')
parser.add_argument('--bns_params', action='store_true', help='Calculate EM parameters using BNS parameters')
args = parser.parse_args()

def convert_time(t0):
    t = Time(t0, format='gps')
    return t.mjd

delta_t = 0
if args.t0 != 0:
    delta_t = 0

if args.time_format == 'gps':
    delta_t = convert_time(delta_t)

### create a dictionary mapping parameter names to values
params = dict(args.p)
for p in params:
    params[p] = float(params[p])

### convert BNS parameters to EM parameters, if needed
### borrowed from EOSManager.py
def lambda_from_m(m):
    eos = lalsim.SimNeutronStarEOSByName("AP4")
    eos_fam = lalsim.CreateSimNeutronStarFamily(eos)
    if m<10**15:
        m=m*lal.MSUN_SI
    k2=lalsim.SimNeutronStarLoveNumberK2(m, eos_fam)
    r=lalsim.SimNeutronStarRadius(m, eos_fam)
    m=m*lal.G_SI/lal.C_SI**2
    lam=2./(3*lal.G_SI)*k2*r**5
    dimensionless_lam=lal.G_SI*lam*(1/m)**5
    return dimensionless_lam

if args.bns_params:
    lambda1 = lambda_from_m(params['m1'])
    lambda2 = lambda_from_m(params['m2'])
    params['mej'] = calc_mej(params['m1'], lambda1, params['m2'], lambda2)
    params['vej'] = calc_vej(params['m1'], lambda1, params['m2'], lambda2)

### initialize the model
if args.orientation is not None:
    profile = args.orientation[0]
    angle = float(args.orientation[1])
    model = model_dict['oriented'](args.m, profile)
    params['angle'] = angle
else:
    model = model_dict[args.m]()

t_bounds = [args.tmin, args.tmax]

### set model parameters
model.set_params(params, t_bounds)

### generate times
tdays = np.linspace(args.tmin, args.tmax, args.n)

data_dict = {}

### generate the data
for band in model.bands:
    m, _ = model.evaluate(tdays, band)
    if 'dist' in params:
        dist = params['dist']
        m += 5*(np.log10(dist*1e6) - 1)
    data = np.empty((4, args.n))
    data[0] = tdays + delta_t
    data[2] = m + np.random.uniform(-1 * args.err, args.err, args.n) # generate errors
    data[3] = np.ones(args.n) * args.err
    data_dict[band] = data

### if not json format, just save as internal format
if not args.json:
    for band in data_dict:
        data = data_dict[band]
        filename = args.out + band + '.txt'
        np.savetxt(filename, data.T)

### json
else:
    photometry = []
    json_dict = {
        'SyntheticEvent': {
            'photometry': photometry
        }
    }
    for band in data_dict:
        data = data_dict[band]
        for i in range(args.n):
            [t, _, m, m_err] = data[:,i]
            photometry.append({
                'band':band,
                'time':t,
                'magnitude':m,
                'e_magnitude':m_err
            })
    with open(args.out + 'SyntheticEvent.json', 'w') as f:
        json.dump(json_dict, f, indent=4)

### save true values

truths = []
for param in model.param_names:
    truths.append(params[param])
filename = args.out + 'test_truths.txt'
np.savetxt(filename, truths)
