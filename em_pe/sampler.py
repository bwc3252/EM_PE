# -*- coding: utf-8 -*-
'''
Generate Posterior Samples
--------------------------
Generate posterior parameter samples according to the given command line
arguments (if used from command line) or function arguments (if called from
another python script)

CLI example::

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
      --ncomp NCOMP    Number of Gaussian components for integrator
'''

from __future__ import print_function
import numpy as np
import argparse
import sys

from models import model_dict, bounds_dict
from integrator_utils import monte_carlo_integrator

try:
    import progressbar
    progress = True
except:
    print('No progressbar')
    progress = False


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
    parser.add_argument('--ncomp', type=int, default=1, help='Number of Gaussian components for integrator')
    return parser.parse_args()

class sampler:
    '''
    Generate posterior samples. This is the function called when using CLI.

    Parameters
    ----------
    data_loc : string
        Directory containing data files
    m : list
        List of [model_name, weight] pairs
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
    '''
    def __init__(self, data_loc, m, files, out, v=True, L_cutoff=0, min_iter=20,
                 max_iter=20, ncomp=1):
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

        ###variables to store
        self.data = None
        self.bands_used = None
        self.models = None
        self.ordered_params = None
        self.bounds = None
        self.t_bounds = None

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

    def _initialize_models(self):
        if self.v:
            print('Initializing models... ', end='')
        ### initialize model objects
        models = [model_dict[name](float(weight)) for [name, weight] in self.m]
        ordered_params = [] # keep track of all parameters used
        bounds = [] # bounds for each parameter
        for model in models:
            for param in model.param_names:
                if param not in ordered_params:
                    ordered_params.append(param)
                    bounds.append(bounds_dict[param])
        t_bounds = [np.inf, -1 * np.inf] # tmin and tmax for each band
        for band in self.bands_used:
            t = self.data[band][0]
            t_bounds[0] = min(min(t), t_bounds[0])
            t_bounds[1] = max(max(t), t_bounds[1])
        if self.v:
            print('finished')
        self.models = models
        self.ordered_params = ordered_params
        self.bounds = bounds
        self.t_bounds = t_bounds

    def _evaluate_lnL(self, params):
        temp_data = {} # used to hold model data and squared model error
        for band in self.bands_used:
            temp_data[band] = [0, 0]
        ### evaluate each model for the params, in the required bands
        for model in self.models:
            model.set_params(params, self.t_bounds)
            for band in model.bands:
                if band in self.bands_used:
                    t = self.data[band][0]
                    m, m_err = model.evaluate(t, band)
                    temp_data[band][0] += m * model.weight
                    temp_data[band][1] += m_err**2
        lnL = 0
        ### calculate lnL
        for band in self.bands_used:
            x = self.data[band][2]
            err = self.data[band][3]
            m = temp_data[band][0]
            m_err_squared = temp_data[band][1]
            diff = x - m
            #print(diff)
            lnL += np.sum(diff**2 / (err**2 + m_err_squared))
        return -0.5 * lnL

    def _integrand(self, samples):
        n, _ = samples.shape
        ret = []
        if self.v and progress:
            bar = progressbar.ProgressBar(maxval=n,
                widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
            bar.start()
        i = 0
        for row in samples:
            if progress and self.v and (i % 10 == 0):
                bar.update(i)
            i += 1
            params = dict(zip(self.ordered_params, row)) # map each parameter's name to its value
            lnL = self._evaluate_lnL(params)
            ret.append(lnL)
        if self.v:
            print()
        ret = np.array(ret).reshape((n, 1))
        ret[np.isnan(ret)] = -1 * np.inf
        return ret

    def _generate_samples(self):
        '''
        Generate posterior samples. This is the function called when using CLI.

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
        if self.v:
            print('Generating posterior samples')
            sys.stdout.flush()
        ### get a tuple of integers as dimensions to make the integrator happy
        param_ind = tuple(range(len(self.ordered_params)))
        dim = len(self.ordered_params) # number of dimensions
        gmm_dict = {param_ind:None}
        ### initialize and run the integrator
        integrator = monte_carlo_integrator.integrator(dim, self.bounds, gmm_dict, self.ncomp,
                        proc_count=None, L_cutoff=self.L_cutoff, use_lnL=True,
                        user_func=sys.stdout.flush())
        integrator.integrate(self._integrand, min_iter=self.min_iter, max_iter=self.max_iter)
        ### make the array of samples
        samples = integrator.cumulative_values
        samples = np.append(samples, integrator.cumulative_p, axis=1)
        samples = np.append(samples, integrator.cumulative_p_s, axis=1)
        samples = np.append(samples, integrator.cumulative_samples, axis=1)
        if self.v:
            print('Integral result:', integrator.integral)
        return samples

    def generate_samples(self):
        '''
        Generate posterior samples.
        '''
        self._read_data()
        self._initialize_models()
        samples = self._generate_samples()
        header = 'lnL p p_s ' + ' '.join(self.ordered_params)
        np.savetxt(self.out, samples, header=header)

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
    s = sampler(data_loc, m, files, out, v, L_cutoff, min_iter, max_iter, ncomp)
    s.generate_samples()

if __name__ == '__main__':
    main()
