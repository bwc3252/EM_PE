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
                          [--truth_file TRUTH_FILE] [--out OUT] [--p P] [--c C]

    Generate corner plot from posterior samples

    optional arguments:
      -h, --help            show this help message and exit
      --posterior_samples POSTERIOR_SAMPLES
                            File with posterior samples
      --truth_file TRUTH_FILE
                            File with true parameter values
      --out OUT             File to save plot to
      --p P                 Parameter name to plot
      --c C                 Minimum likelihood for points to keep
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
    parser.add_argument('--c', type=float, default=0, help='Minimum likelihood for points to keep. Takes precedence over --frac')
    parser.add_argument('--frac', type=float, default=1.0, help='Fraction of points to keep')
    parser.add_argument('--legend', action='append', help='Name of posterior sample set for plot legend. Assumed to be in the same order as the posterior sample files')
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
    if args.c <= 0:
        min_lnL = -1 * np.inf
    else:
        min_lnL = np.log(args.c)
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
        lnL = samples[:,0]
        p = samples[:,1]
        p_s = samples[:,2]
        if args.c != 0:
            mask = lnL > min_lnL
        elif args.frac != 1.0:
            ind = np.argsort(lnL)
            n = int(len(lnL) * args.frac)
            mask = ind[len(lnL) - n:]
        else:
            mask = [True] * len(lnL)
        lnL = lnL[mask]
        p = p[mask]
        p_s = p_s[mask]
        ### get columns of array corresponding to actual parameter samples
        x = samples[:,[index_dict[name] for name in args.p]]
        lnL += abs(np.max(lnL))
        L = np.exp(lnL)
        weights = L * p / p_s
        weights /= np.sum(weights)
        mask2 = weights > 0
        print(np.sum(mask2), 'samples with weight > 0')
        weights = weights[mask2]
        x = x[mask]
        x = x[mask2]
        color = color_list[i % len(color_list)]
        #levels = None
        levels = [0.5, 0.9]
        fig_base = corner.corner(x, weights=weights, levels=levels, fig=fig_base, labels=param_names, truths=truths,
                                 color=color, plot_datapoints=False, plot_density=False, no_fill_contours=True,
                                 contours=True)
        i += 1
    if args.legend != []:
        xcoord = len(args.p)
        ycoord = len(args.p)
        lgd = plt.legend(args.legend, bbox_to_anchor=(xcoord, ycoord), loc="center right")
    plt.savefig(args.out, bbox_extra_artists=(lgd,), bbox_inches='tight')

if __name__ == '__main__':
    generate_plot()
