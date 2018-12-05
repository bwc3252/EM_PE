# -*- coding: utf-8 -*-
'''
Plot corner
-----------
Generates and saves a corner plot from posterior samples.

Example::

    $ python plot_corner.py --posterior_samples samples.txt --truth_file truths.txt --out fig.png --p a --p b

To see full command line parameter documentation::

    $ python plot_corner.py -h
    usage: plot_corner.py [-h] [--posterior_samples POSTERIOR_SAMPLES]
                          [--truth_file TRUTH_FILE] [--out OUT] [--p P]

    Generate corner plot from posterior samples

    optional arguments:
      -h, --help            show this help message and exit
      --posterior_samples POSTERIOR_SAMPLES
                            File with posterior samples
      --truth_file TRUTH_FILE
                            File with true parameter values
      --out OUT             File to save plot to
      --p P                 Parameter name to plot
'''

from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import corner
import argparse

def _parse_command_line_args():
    '''
    Parses and returns the command line arguments.
    '''
    parser = argparse.ArgumentParser(description='Generate corner plot from posterior samples')
    parser.add_argument('--posterior_samples', action='append', help='File with posterior samples')
    parser.add_argument('--truth_file', help='File with true parameter values')
    parser.add_argument('--out', help='File to save plot to')
    parser.add_argument('--p', action='append', help='Parameter name to plot')
    return parser.parse_args()

def generate_plot():
    '''
    Generates a corner plot for the specified posterior samples and parameters.
    '''
    args = _parse_command_line_args()
    ### colors to iterate through
    color_list=['black', 'red', 'green', 'blue','yellow']
    sample_files = args.posterior_samples
    truth_file = args.truth_file
    if truth_file is not None:
        truths = np.loadtxt(truth_file)
    else:
        truths = None
    fig_base = None
    i = 0
    for file in sample_files:
        samples = np.loadtxt(file, skiprows=1)
        with open(file) as f:
            header = f.readline().strip().split(' ')
        param_names = header[4:]
        index_dict = {}
        for index in range(4, len(header)):
            index_dict[header[index]] = index - 1
        L = samples[:,0]
        p = samples[:,1]
        p_s = samples[:,2]
        n, m = samples.shape
        ### get columns of array corresponding to actual parameter samples
        x = samples[:,[index_dict[name] for name in args.p]]
        weights = L * p / p_s
        weights /= np.sum(weights)
        color = color_list[i % len(color_list)]
        fig_base = corner.corner(x, weights=weights, fig=fig_base, labels=param_names, truths=truths,
                                 color=color, plot_datapoints=False, plot_density=False, no_fill_contours=True,
                                 contours=True,)
        i += 1
    plt.savefig(args.out)

if __name__ == '__main__':
    generate_plot()
