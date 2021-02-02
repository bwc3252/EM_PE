# -*- coding: utf-8 -*-

from __future__ import print_function
import numpy as np
import argparse
import sys
from multiprocessing import Pool

### hacky fix because the import system is different when run as a package vs.
### run as a script
try:
    from models import model_dict, param_dict
except ModuleNotFoundError:
    from .models import model_dict, param_dict

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
    parser.add_argument('--ncomp', action='append', nargs=2, help='Number of Gaussian components for a given dimension')
    parser.add_argument('--fixed-param', action='append', nargs=2, help='Parameters with fixed values')
    parser.add_argument('--estimate-dist', action="store_true", help='Estimate distance')
    parser.add_argument('--epoch', type=int, default=100, help='Iterations before resetting sampling distributions')
    parser.add_argument('--correlate-dims', action='append', nargs='+', help='Parameters to group together')
    parser.add_argument('--burn-in', type=int, help='Number or iterations to use as burn-in')
    parser.add_argument('--beta-start', type=float, default=1.0, help='Initial beta value: if burn-in is set, 0 < beta <= 1 is multiplied by lnL for every iteration of burn-in')
    parser.add_argument('--beta-end', type=float, default=1.0, help='')
    parser.add_argument('--keep-npts', type=int, help='Store the n highest-likelihood samples')
    parser.add_argument('--nprocs', type=int, default=1, help='Number of parallel processes to use for likelihood evaluation')
    parser.add_argument('--set-limit', action='append', nargs=3, help='Modify parameter limits (e.g. --set-limit mej_red 0.008 0.012)')
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
                 max_iter=20, ncomp=None, fixed_params=None,
                 estimate_dist=True, epoch=5, correlate_dims=None, burn_in_length=None,
                 beta_start=1.0, beta_end=1.0, keep_npts=None, nprocs=1, limits=None):
        ### parameters passed in from user or main()
        self.data_loc = data_loc
        self.m = m
        self.files = files
        self.out = out
        self.v = v
        self.L_cutoff = L_cutoff
        self.min_iter = min_iter
        self.max_iter = max_iter
        self.estimate_dist = estimate_dist
        self.epoch = epoch
        self.fixed_params = {}
        self.correlate_dims = correlate_dims
        self.burn_in_length = burn_in_length
        #self.burn_in_start = burn_in_start
        #self.burn_in_npts = burn_in_start
        self.beta_start = beta_start
        self.beta_end = beta_end
        self.keep_npts = keep_npts
        self.nprocs = nprocs
        self.limits = limits if limits is not None else {}
        if ncomp is None:
            self.ncomp = 1
        else:
            self.ncomp = {name:int(nc) for [name, nc] in ncomp}

        ### convert types for fixed params, make it a dict
        if fixed_params is not None:
            for [name, value] in fixed_params:
                self.fixed_params[name] = float(value)

        ###variables to store
        self.integrator = None
        self.data = None
        self.bands_used = None
        self.models = []
        self.params = None
        self.ordered_params = None
        self.bounds = None
        self.t_bounds = None
        self.iteration = 0
        self.iteration_size = 0

        self.cumulative_lnL = np.array([])

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
            params = {}
            for param in model.param_names:
                if param not in ordered_params and param not in self.fixed_params:
                    ordered_params.append(param)
                    params[param] = param_dict[param]()
                    if param in self.limits.keys():
                        llim, rlim = self.limits[param]
                        params[param].update_limits(llim, rlim)
                    bounds.append([params[param].llim, params[param].rlim])
        t_bounds = [np.inf, -1 * np.inf] # tmin and tmax for each band
        for band in self.bands_used:
            t = self.data[band][0]
            t_bounds[0] = min(min(t), t_bounds[0])
            t_bounds[1] = max(max(t), t_bounds[1])
        if self.estimate_dist:
            ordered_params.append('dist')
            params["dist"] = param_dict["dist"]()
            bounds.append([params["dist"].llim, params["dist"].rlim])
        self.params = params
        self.ordered_params = ordered_params
        self.bounds = bounds
        self.t_bounds = t_bounds
        if self.v:
            print('finished')

    def _prior(self, sample_array):
        n, m = sample_array.shape
        ret = np.ones(n)
        index_dict = dict(zip(self.ordered_params, range(len(self.ordered_params))))
        ### evaluate the prior for each parameter
        for p in self.ordered_params:
            x = sample_array[:,index_dict[p]]
            ret *= self.params[p].prior(x)
        return ret.reshape((n, 1))

    def _evaluate_lnL(self, params, model, vectorized=False):
        temp_data = {} # used to hold model data and squared model error
        for band in self.bands_used:
            temp_data[band] = [None, None]
        ### evaluate each model for the params, in the required bands
        model.set_params(params, self.t_bounds)
        for band in model.bands:
            if band not in self.bands_used:
                continue
            t = self.data[band][:,0]
            m, m_err = model.evaluate(t, band)
            temp_data[band][0] = m
            temp_data[band][1] = m_err
            if 'dist' in params:
                dist = params['dist']
                if vectorized:
                    for i in range(dist.size):
                        temp_data[band][0][i] += 5.0 * (np.log10(dist[i] * 1.0e6) - 1.0)
                else:
                    temp_data[band][0] += 5.0 * (np.log10(dist * 1.0e6) - 1.0)

        if vectorized:
            for band in self.bands_used:
                x = self.data[band][:,2]
                err = self.data[band][:,3]
                m = temp_data[band][0]
                m_err = temp_data[band][1]
                lnL = np.empty(m.shape[0])
                for i in range(m.shape[0]):
                    diff = x - m[i]
                    lnL[i] += np.sum(diff**2 / (err**2 + m_err[i]**2) + np.log(2.0 * np.pi * (err**2 + m_err[i]**2)))
            return -0.5 * lnL

        lnL = 0
        ### calculate lnL
        for band in self.bands_used:
            x = self.data[band][:,2]
            err = self.data[band][:,3]
            m = temp_data[band][0]
            m_err = temp_data[band][1]
            diff = x - m
            lnL += np.sum(diff**2 / (err**2 + m_err**2) + np.log(2.0 * np.pi * (err**2 + m_err**2)))
        return -0.5 * lnL

    def _get_current_samples(self):
        samples = self.cumulative_lnL.reshape((self.cumulative_lnL.size, 1))
        samples = np.append(samples, self.integrator.cumulative_p, axis=1)
        samples = np.append(samples, self.integrator.cumulative_p_s, axis=1)
        samples = np.append(samples, self.integrator.cumulative_samples, axis=1)
        if self.keep_npts is not None and self.keep_npts < samples.shape[0]:
            ind_sorted = np.argsort(samples[:,0])
            samples = samples[ind_sorted[samples.shape[0] - self.keep_npts:]]
        return samples

    def _integrand_subprocess(self, arg):
        model, samples = arg

        ### if the model is vectorized, use that
        if model.vectorized:
            params = {}
            for i, p in enumerate(self.ordered_params):
                params[p] = samples[:,i]
            for p in self.fixed_params:
                params[p] = self.fixed_params[p] * np.ones(samples.shape[0])
            return self._evaluate_lnL(params, model, vectorized=True)

        ### otherwise do it in a loop
        ret = np.empty(samples.shape[0])
        for i in range(samples.shape[0]):
            row = samples[i]
            params = dict(zip(self.ordered_params, row)) # map each parameter's name to its value
            for p in self.fixed_params:
                params[p] = self.fixed_params[p]
            ret[i] = self._evaluate_lnL(params, model)
        return ret

    def _integrand(self, samples):
        if self.v:
            print("Iteration", self.iteration)
        #if self.iteration > 0:
        #    header = 'lnL p p_s ' + ' '.join(self.ordered_params)
        #    if self.v:
        #        print("saving intermediate samples...")
        #    fname = self.out.split(".")[0] + "_intermediate." + "".join(self.out.split(".")[1:])
        #    np.savetxt(fname, self._get_current_samples(), header=header)
        if self.burn_in_length is not None and self.iteration < self.burn_in_length:
            beta = np.exp((1.0 - self.iteration / (self.burn_in_length + 1.0)) * np.log(self.beta_start)
                    + self.iteration * np.log(self.beta_end) / (self.burn_in_length + 1.0)) # evenly-spaced on log scale
        else:
            beta = self.beta_end
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
        #print(np.min(ret), np.max(ret))
        ret[np.isnan(ret)] = -1 * np.inf
        self.iteration += 1
        self.cumulative_lnL = np.append(self.cumulative_lnL, ret)
        ret *= beta
        if self.v:
            print("points with non-zero likelihood:", np.sum(np.exp(ret - np.max(ret)) > 0.0))
        return ret

    def _generate_samples(self):
        if self.v:
            print('Generating posterior samples')
            sys.stdout.flush()
        dim = len(self.ordered_params) # number of dimensions
        ### dictionary mapping parameter names to ordered_params index
        param_ind = dict(zip(self.ordered_params, range(dim)))
        if self.correlate_dims is None:
            ### assume separately-sampled dimensions
            gmm_dict = {(i,):None for i in range(dim)}
        else:
            ### make a flattened version of correlate_dims
            correlated_params = [item for sublist in self.correlate_dims for item in sublist]
            ### add the correlated dimensions to the gmm_dict
            gmm_dict = {tuple([param_ind[p] for p in sublist]):None for sublist in self.correlate_dims}
            ### add the non-correlated dimensions
            for p in self.ordered_params:
                if p not in correlated_params:
                    gmm_dict[(param_ind[p],)] = None
        if self.ncomp == 1:
            ncomp = 1
        else: # take care of multi-component GMMs
            ### basically just copy gmm_dict (except map tuples to 1 instead of None)
            ncomp = {}
            for ind_tuple in gmm_dict:
                ncomp[ind_tuple] = 1
            for p in self.ncomp.keys():
                ind = param_ind[p]
                for ind_tuple in ncomp.keys():
                    if ind in ind_tuple:
                        ncomp[ind_tuple] = self.ncomp[p]
                        break
        ### initialize and run the integrator
        self.integrator = monte_carlo_integrator.integrator(dim, self.bounds, gmm_dict, ncomp,
                        proc_count=None, L_cutoff=self.L_cutoff, use_lnL=True,
                        user_func=sys.stdout.flush(), prior=self._prior)
        self.integrator.integrate(self._integrand, min_iter=self.min_iter, max_iter=self.max_iter, 
                progress=self.v, epoch=self.epoch)
        ### make the array of samples
        if self.v:
            print('Integral result:', self.integrator.integral)
        return self._get_current_samples()

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
    #if args.burn_in is not None:
    #    burn_in_length = int(args.burn_in[0])
    #    burn_in_start = int(args.burn_in[1])
    #else:
    #    burn_in_length = None
    #    burn_in_start = None
    burn_in_length = args.burn_in
    beta_start = args.beta_start
    beta_end = args.beta_end
    keep_npts = args.keep_npts
    nprocs = args.nprocs
    if args.set_limit is not None:
        limits = {name:(float(llim), float(rlim)) for (name, llim, rlim) in args.set_limit}
    else:
        limits = None
    s = sampler(data_loc, m, files, out, v, L_cutoff, min_iter, max_iter, ncomp, 
            fixed_params, estimate_dist, epoch, correlate_dims,
            burn_in_length, beta_start, beta_end, keep_npts, nprocs, limits)
    #        burn_in_length, burn_in_start, beta_start, keep_npts, nprocs)
    s.generate_samples()

if __name__ == '__main__':
    main()
