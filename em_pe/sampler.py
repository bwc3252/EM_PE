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
    parser.add_argument('--m', help='Name of model to use')
    parser.add_argument('-v', action='store_true', help='Verbose mode')
    parser.add_argument('--cutoff', default=0, type=float, help='Likelihood cutoff for storing posterior samples')
    parser.add_argument('--f', action='append', help='Name of a data file')
    parser.add_argument('--min', default=20, type=int, help='Minimum number of integrator iterations')
    parser.add_argument('--max', default=20, type=int, help='Maximum number of integrator iterations')
    parser.add_argument('--out', help='Location to store posterior samples')
    parser.add_argument('--ncomp', type=int, default=1, help='Number of Gaussian components for integrator')
    parser.add_argument('--fixed_param', action='append', nargs=2, help='Parameters with fixed values')
    parser.add_argument('--orientation', help='Orientation dependance to use (defaults to None)')
    parser.add_argument('--estimate_dist', action="store_true", help='Estimate distance')
    parser.add_argument('--epoch', type=int, default=100, help='Iterations before resetting sampling distributions')
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
                 max_iter=20, ncomp=1, fixed_params=None, orientation=None,
                 estimate_dist=True, epoch=5):
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
        self.orientation = orientation
        self.estimate_dist = estimate_dist
        self.epoch = epoch
        self.fixed_params = {}

        ### convert types for fixed params, make it a dict
        if fixed_params is not None:
            for [name, value]  in fixed_params:
                self.fixed_params[name] = float(value)

        ###variables to store
        self.data = None
        self.bands_used = None
        self.models = None
        self.ordered_params = None
        self.bounds = None
        self.t_bounds = None

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
        ### initialize model object
        if self.orientation is not None:
            model = model_dict['oriented'](self.m, self.orientation)
        else:
            model = model_dict[self.m]()
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
        self.model = model
        self.ordered_params = ordered_params
        self.bounds = bounds
        self.t_bounds = t_bounds

    def _evaluate_lnL(self, params):
        temp_data = {} # used to hold model data and squared model error
        for band in self.bands_used:
            temp_data[band] = [0, 0]
        ### evaluate each model for the params, in the required bands
        self.model.set_params(params, self.t_bounds)
        for band in self.model.bands:
            if band in self.bands_used:
                t = self.data[band][:,0]
                m, m_err = self.model.evaluate(t, band)
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
            ### take care of fixed parameters
            for p in self.fixed_params:
                params[p] = self.fixed_params[p]
            lnL = self._evaluate_lnL(params)
            ret.append(lnL)
        if self.v:
            print()
        ret = np.array(ret).reshape((n, 1))
        ret[np.isnan(ret)] = -1 * np.inf
        return ret

    def _generate_samples(self):
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
        integrator.integrate(self._integrand, min_iter=self.min_iter, max_iter=self.max_iter, 
                progress=self.v, epoch=self.epoch)
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
        lnL = np.empty(n)
        row = 0
        while row < n:
            params = dict(zip(self.ordered_params, sample_array[row]))
            for p in self.fixed_params:
                params[p] = self.fixed_params[p]
            lnL[row] = self._evaluate_lnL(params)
            row += 1
        return lnL

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
    orientation = args.orientation
    estimate_dist = args.estimate_dist
    epoch = args.epoch
    s = sampler(data_loc, m, files, out, v, L_cutoff, min_iter, max_iter, ncomp, 
            fixed_params, orientation, estimate_dist, epoch)
    s.generate_samples()

if __name__ == '__main__':
    main()
