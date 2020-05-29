# -*- coding: utf-8 -*-
'''
Generate Posterior Samples
--------------------------
Generate posterior parameter samples according to the given command line
arguments (if used from command line) or function arguments (if called from
another python script)

CLI Usage
^^^^^^^^^
Example::

    $ python sampler.py --dat [data directory] --m [model name] --f [data filename] --out [file to save samples to]

To see full command line parameter documentation::

    $ python sampler.py -h
    usage: sampler.py [-h] [--dat DAT] [--m M] [-v] [--cutoff CUTOFF] [--f F]
                      [--min MIN] [--max MAX] [--out OUT] [--ncomp NCOMP]
                      [--fixed_param FIXED_PARAM FIXED_PARAM]

    Generate posterior parameter samples from lightcurve data

    optional arguments:
      -h, --help            show this help message and exit
      --dat DAT             Path to data directory
      --m M                 Name of model to use
      -v                    Verbose mode
      --cutoff CUTOFF       Likelihood cutoff for storing posterior samples
      --f F                 Name of a data file
      --min MIN             Minimum number of integrator iterations
      --max MAX             Maximum number of integrator iterations
      --out OUT             Location to store posterior samples
      --ncomp NCOMP         Number of Gaussian components for integrator
      --fixed_param FIXED_PARAM FIXED_PARAM
                            Parameters with fixed values

Python API
^^^^^^^^^^
Example::

    from em_pe import sampler

    dat = './'
    m = 'two_comp'
    f = ['H.txt', 'J.txt']
    out = 'posterior_samples.txt'

    ### Initialize sampler
    s = sampler(dat, m, f, out)

    ### Generate samples
    s.generate_samples()
'''

from __future__ import print_function
import numpy as np
import argparse
import sys
from multiprocessing import Pool

### hacky fix because the import system is different when run as a package vs.
### run as a script
try:
    from models import model_dict, bounds_dict
except ModuleNotFoundError:
    from .models import model_dict, bounds_dict

import RIFT.integrators.MonteCarloEnsemble as monte_carlo_integrator

def _parse_command_line_args():
    '''
    Parses and returns the command line arguments.
    '''
    parser = argparse.ArgumentParser(description='Generate posterior parameter samples from lightcurve data')
    parser.add_argument('--dat', help='Path to data directory')
    parser.add_argument('--m', help='Name of model to use')
    parser.add_argument('-v', action='store_true', help='Verbose mode')
    parser.add_argument('--cutoff', default=0, type=float, help='Likelihood cutoff for storing posterior samples')
    parser.add_argument('--f', action='append', help='Name of a data file')
    parser.add_argument('--min', default=20, type=int, help='Minimum number of integrator iterations')
    parser.add_argument('--max', default=20, type=int, help='Maximum number of integrator iterations')
    parser.add_argument('--out', help='Location to store posterior samples')
    parser.add_argument('--ncomp', type=int, default=1, help='Number of Gaussian components for integrator')
    parser.add_argument('--fixed-param', action='append', nargs=2, help='Parameters with fixed values')
    parser.add_argument('--estimate-dist', action="store_true", help='Estimate distance')
    parser.add_argument('--epoch', type=int, default=100, help='Iterations before resetting sampling distributions')
    parser.add_argument('--correlate-dims', action='append', nargs='+', help='Parameters to group together')
    parser.add_argument('--burn-in', nargs=2, help='Number or iterations to use as burn-in and number of samples to start with per band')
    parser.add_argument('--beta-start', type=float, default=1.0, help='Initial beta value: if burn-in is set, 0 < beta <= 1 is multiplied by lnL for every iteration of burn-in')
    parser.add_argument('--keep-npts', type=int, help='Store the n highest-likelihood samples')
    parser.add_argument('--nprocs', type=int, default=1, help='Number of parallel processes to use for likelihood evaluation')
    return parser.parse_args()

class sampler:
    '''
    Generate posterior samples. This is the function called when using CLI.

    Parameters
    ----------
    data_loc : string
        Directory containing data files
    m : string
        Name of model to use
    files : list
        List of data file names (in the data_loc directory)
    v : bool
        Verbose
    L_cutoff : float
        Likelihood cutoff for saving samples
    min_iter : int
        Minimum number of integrator iterations
    max_iter : int
        Maximum number of integrator iterations
    ncomp : int
        Number of Gaussian components to use for integrator
    fixed_params : list
        List of [param_name, value] pairs
    '''
    def __init__(self, data_loc, m, files, out, v=True, L_cutoff=0, min_iter=20,
                 max_iter=20, ncomp=1, fixed_params=None,
                 estimate_dist=True, epoch=5, correlate_dims=None, burn_in_length=None,
                 burn_in_start=None, beta_start=1.0, keep_npts=None, nprocs=1):
        ### parameters passed in from user or main()
        self.data_loc = data_loc
        self.m = m
        self.files = files
        self.out = out
        self.v = v
        self.L_cutoff = L_cutoff
        self.min_iter = min_iter
        self.max_iter = max_iter
        self.ncomp = ncomp
        self.estimate_dist = estimate_dist
        self.epoch = epoch
        self.fixed_params = {}
        self.correlate_dims = correlate_dims
        self.burn_in_length = burn_in_length
        self.burn_in_start = burn_in_start
        self.burn_in_npts = burn_in_start
        self.beta_start = beta_start
        self.keep_npts = keep_npts
        self.nprocs = nprocs

        ### convert types for fixed params, make it a dict
        if fixed_params is not None:
            for [name, value]  in fixed_params:
                self.fixed_params[name] = float(value)

        ###variables to store
        self.data = None
        self.bands_used = None
        self.models = []
        self.ordered_params = None
        self.bounds = None
        self.t_bounds = None
        self.iteration = 0
        self.iteration_size = 0

        ### initialization things
        self._read_data()
        self._initialize_model()

    def _read_data(self):
        if self.v:
            print('Loading data... ', end='')
        data = {}
        bands_used = []
        for fname in self.files:
            band = fname.split('.txt')[0]
            data[band] = np.loadtxt(self.data_loc + fname)
            bands_used.append(band)
        if self.v:
            print('finished')
        self.data = data
        self.bands_used = bands_used

    def _initialize_model(self):
        if self.v:
            print('Initializing models... ', end='')
        ### initialize model objects (one for each parallel process)
        for i in range(self.nprocs):
            model = model_dict[self.m]()
            self.models.append(model)
            ordered_params = [] # keep track of all parameters used
            bounds = [] # bounds for each parameter
            for param in model.param_names:
                if param not in ordered_params and param not in self.fixed_params:
                    ordered_params.append(param)
                    bounds.append(bounds_dict[param])
        t_bounds = [np.inf, -1 * np.inf] # tmin and tmax for each band
        for band in self.bands_used:
            t = self.data[band][0]
            t_bounds[0] = min(min(t), t_bounds[0])
            t_bounds[1] = max(max(t), t_bounds[1])
        if self.estimate_dist:
            ordered_params.append('dist')
            bounds.append(bounds_dict['dist'])
        if self.v:
            print('finished')
        self.ordered_params = ordered_params
        self.bounds = bounds
        self.t_bounds = t_bounds

    def _evaluate_lnL(self, params, model):
        temp_data = {} # used to hold model data and squared model error
        for band in self.bands_used:
            temp_data[band] = [0, 0]
        ### evaluate each model for the params, in the required bands
        model.set_params(params, self.t_bounds)
        for band in model.bands:
            if band in self.bands_used:
                t = self.data[band][:,0]
                m, m_err = model.evaluate(t, band)
                temp_data[band][0] += m
                temp_data[band][1] += m_err**2
                if 'dist' in params:
                    dist = params['dist']
                    temp_data[band][0] += 5*(np.log10(dist*1e6) - 1)
        lnL = 0
        ### calculate lnL
        for band in self.bands_used:
            x = self.data[band][:,2]
            err = self.data[band][:,3]
            m = temp_data[band][0]
            m_err_squared = temp_data[band][1]
            diff = x - m
            if self.burn_in_length is not None and self.iteration < self.burn_in_length:
                npts = int(self.burn_in_start + (self.iteration / self.burn_in_length) * (diff.size - self.burn_in_start))
                if npts < diff.size:
                    ind = np.random.choice(diff.size, size=self.burn_in_npts, replace=False)
                    diff = diff[ind]
                    err = err[ind]
                    #m_err_squared = m_err_squared[ind] # FIXME fix this when the models start actually having error bars
            lnL += np.sum(diff**2 / (err**2 + m_err_squared))
        return -0.5 * lnL

    def _integrand_subprocess(self, arg):
        model, samples = arg
        ret = np.empty(samples.shape[0])
        for i in range(samples.shape[0]):
            row = samples[i]
            params = dict(zip(self.ordered_params, row)) # map each parameter's name to its value
            for p in self.fixed_params:
                params[p] = self.fixed_params[p]
            ret[i] = self._evaluate_lnL(params, model)
        return ret

    def _integrand(self, samples):
        beta = self.beta_start + (1.0 - self.beta_start) * (self.iteration / self.burn_in_length)
        n, _ = samples.shape
        self.iteration_size = n
        if self.nprocs == 1:
            ret = self._integrand_subprocess((self.models[0], samples))
        elif self.nprocs > 1:
            samples_split = np.array_split(samples, self.nprocs)
            args = zip(self.models, samples_split)
            with Pool(self.nprocs) as p:
                ret = np.concatenate(p.map(self._integrand_subprocess, args))
        else:
            raise RuntimeError("nprocs < 1: How can you have less than 1 process?")
        ret = ret.reshape((n, 1))
        ret[np.isnan(ret)] = -1 * np.inf
        self.iteration += 1
        return ret * beta

    def _generate_samples(self):
        if self.v:
            print('Generating posterior samples')
            sys.stdout.flush()
        dim = len(self.ordered_params) # number of dimensions
        if self.correlate_dims is None:
            ### assume separately-sampled dimensions
            gmm_dict = {(i,):None for i in range(dim)}
        else:
            ### dictionary mapping parameter names to ordered_params index
            param_ind = dict(zip(self.ordered_params, range(dim)))
            ### make a flattened version of correlate_dims
            correlated_params = [item for sublist in self.correlate_dims for item in sublist]
            ### add the correlated dimensions to the gmm_dict
            gmm_dict = {tuple([param_ind[p] for p in sublist]):None for sublist in self.correlate_dims}
            ### add the non-correlated dimensions
            for p in self.ordered_params:
                if p not in correlated_params:
                    gmm_dict[(param_ind[p],)] = None
        ### initialize and run the integrator
        integrator = monte_carlo_integrator.integrator(dim, self.bounds, gmm_dict, self.ncomp,
                        proc_count=None, L_cutoff=self.L_cutoff, use_lnL=True,
                        user_func=sys.stdout.flush())
        integrator.integrate(self._integrand, min_iter=self.min_iter, max_iter=self.max_iter, 
                progress=self.v, epoch=self.epoch)
        ### make the array of samples
        samples = integrator.cumulative_values
        samples = np.append(samples, integrator.cumulative_p, axis=1)
        samples = np.append(samples, integrator.cumulative_p_s, axis=1)
        samples = np.append(samples, integrator.cumulative_samples, axis=1)
        if self.burn_in_length is not None:
            samples = samples[self.burn_in_length * self.iteration_size:]
        if self.keep_npts is not None and self.keep_npts < samples.shape[0]:
            ind_sorted = np.argsort(samples[:,0])
            samples = samples[ind_sorted[samples.shape[0] - self.keep_npts:]]
        if self.v:
            print('Integral result:', integrator.integral)
        return samples

    def generate_samples(self):
        '''
        Generate posterior samples.
        '''
        samples = self._generate_samples()
        header = 'lnL p p_s ' + ' '.join(self.ordered_params)
        np.savetxt(self.out, samples, header=header)

    def log_likelihood(self, samples, vect=False):
        '''
        Function to calculate log-likelihoods of parameter samples

        Parameters
        ----------
        samples : dict
            Dictionary mapping parameter names to values. The values can either
            be single floats or arrays. If arrays, the vect parameter should be
            True
        vect : bool
            Set to True if passing arrays of parameter samples

        Returns
        -------
        np.ndarray
            Array of log-likelihoods
        '''
        ### get samples into common format
        for param in samples:
            if vect:
                samples[param] = samples[param].flatten()
                n = samples[param].size
            else:
                samples[param] = np.array(samples[param])
                n = 1
        ### concatenate everything into a single array (so we can take slices of it)
        sample_array = np.empty((n, len(self.ordered_params)))
        for col in range(len(self.ordered_params)):
            sample_array[:,col] = samples[self.ordered_params[col]]
        return self._integrand(sample_array)

def main():
    args = _parse_command_line_args()
    data_loc = args.dat
    m = args.m
    v = args.v
    L_cutoff = args.cutoff
    files = args.f
    min_iter = args.min
    max_iter = args.max
    out = args.out
    ncomp = args.ncomp
    fixed_params = args.fixed_param
    estimate_dist = args.estimate_dist
    epoch = args.epoch
    correlate_dims = args.correlate_dims
    if args.burn_in is not None:
        burn_in_length = int(args.burn_in[0])
        burn_in_start = int(args.burn_in[1])
    beta_start = args.beta_start
    keep_npts = args.keep_npts
    nprocs = args.nprocs
    s = sampler(data_loc, m, files, out, v, L_cutoff, min_iter, max_iter, ncomp, 
            fixed_params, estimate_dist, epoch, correlate_dims,
            burn_in_length, burn_in_start, beta_start, keep_npts, nprocs)
    s.generate_samples()

if __name__ == '__main__':
    main()
