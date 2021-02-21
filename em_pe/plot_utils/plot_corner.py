# -*- coding: utf-8 -*-

from __future__ import print_function
import numpy as np
import corner
import argparse

try:
    import matplotlib.pyplot as plt
except:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

def _parse_command_line_args():
    '''
    Parses and returns the command line arguments.
    '''
    parser = argparse.ArgumentParser(description='Generate corner plot from posterior samples')
    parser.add_argument('--posterior-samples', action='append', help='File with posterior samples')
    parser.add_argument('--truth-file', help='File with true parameter values')
    parser.add_argument('--out', help='File to save plot to')
    parser.add_argument('--p', action='append', help='Parameter name to plot')
    parser.add_argument('--c', type=float, default=0, help='Minimum likelihood for points to keep. Takes precedence over --frac')
    parser.add_argument('--frac', type=float, default=1.0, help='Fraction of points to keep')
    parser.add_argument('--legend', action='append', help='Name of posterior sample set for plot legend. Assumed to be in the same order as the posterior sample files')
    parser.add_argument('--title', help='Title for plot')
    parser.add_argument('--cl', action='append', type=float, help='Adds a confidence level. Set to "default" for default contours')
    parser.add_argument('--combine', action='store_true', help='Generate a plot using ALL sample files specified')
    parser.add_argument('--min-weight', type=float, default=0.0, help='Minimum weight to keep')
    return parser.parse_args()

def generate_corner_plot(sample_files, out, params, truths=None, cutoff=0, frac=1.0, leg=None,
                  cl='default', title=None, combine=False, min_weight=0):
    '''
    Generates a corner plot for the specified posterior samples and parameters.

    Parameters
    ----------
    sample_files : list
        List of posterior sample files
    out : string
        File to save plot to
    params : list
        List of parameter names
    truths : string
        File with true parameter values
    cutoff : float
        Minimum likelihood for points to keep. Takes precedence over frac
    frac : float
        Fraction of points to keep
    leg : list
        List of names of posterior sample set for plot legend. Assumed to be in
        the same order as the posterior sample files
    cl : string or list
        List of confidence levels to plot. Set to 'default' for default contours.
    title : str
        Title for plot
    combine : bool
        Combine samples for all bands into a single dataset and plot this in
        addition to individual bands (useful for generating density gradient)
    '''
    ### colors to iterate through
    color_list=['black', 'red', 'orange', 'yellow', 'green', 'cyan', 'blue',
                'purple', 'gray']
    ### dictionary of LaTeX parameter name strings
    tex_dict = {'mej':'$m_{ej}$ $(M_\\odot)$',
                #'log_mej':'$\\log(m_{ej})$',
                #'log_mej_blue':'$\\log(m_{ej} / M_\\odot)$ (blue)',
                #'log_mej_purple':'$\\log(m_{ej} / M_\\odot)$ (purple)',
                #'log_mej_red':'$\\log(m_{ej} / M_\\odot)$ (red)',
                'mej1':'$m_{ej1}$ $(M_\\odot)$',
                'mej2':'$m_{ej2}$ $(M_\\odot)$',
                'vej':'$v_{ej} / c$',
                'vej1':'$v_{ej1}$ $(v/c)$',
                'vej2':'$v_{ej2}$ $(v/c)$',
                'mej_red':'$m_{ej}$ [red] $(M_\\odot)$',
                'mej_purple':'$m_{ej}$ [purple] $(M_\\odot)$',
                'mej_blue':'$m_{ej}$ [blue] $(M_\\odot)$',
                'mej_dyn':'$m_{ej}$ [dyn] $(M_\\odot)$',
                'mej_wind':'$m_{ej}$ [wind] $(M_\\odot)$',
                'vej_red':'$v_{ej} / c$ [red]',
                'vej_purple':'$v_{ej} / c$ [purple]',
                'vej_blue':'$v_{ej} / c$ [blue]',
                'vej_dyn':'$v_{ej} / c$ [dyn]',
                'vej_wind':'$v_{ej} / c$ [wind]',
                'Tc_red':'Tc [red]',
                'Tc_purple':'Tc [purple]',
                'Tc_blue':'Tc [blue]',
                'm1':'$m_1$ $(M_\\odot)$',
                'm2':'$m_2$ $(M_\\odot)$',
                'dist':'Distance (Mpc)',
                'sigma':'$\\sigma$ (mag)',
                'kappa':'$\kappa$ (cm$^2$ g$^{-1}$)',
                'theta':'$\\theta$ (deg)'
               }
    if cutoff <= 0:
        min_lnL = -1 * np.inf
    else:
        min_lnL = np.log(cutoff)
    if truths is not None:
        with open(truths) as f:
            truths_header = f.readline().strip().split(' ')[1:]
        truths_full = np.loadtxt(truths, skiprows=1)
        truths_dict = {truths_header[i]:truths_full[i] for i in range(truths_full.size)}
    else:
        truths_dict = None
    fig_base = None
    i = 0
    total_samples = []
    headers = []
    for file in sample_files:
        samples = np.loadtxt(file, skiprows=1)
        total_samples.append(samples)
        with open(file) as f:
            ### the "header" contains the column names
            header = f.readline().strip().split(' ')
        headers.append(header)
    if combine:
        total_samples.append(np.concatenate(total_samples, axis=0))
        headers.append(headers[0]) # is this safe?
        if leg is not None:
            leg.append('combined')
    for ind in range(len(total_samples)):
        samples = total_samples[ind]
        header = headers[ind]
        ### the parameter samples are in columns 4 and up, so to get their
        ### names look at the corresponding words in the header
        param_names = header[4:]
        if truths_dict is not None:
            truths = np.array([truths_dict[name] for name in params])
        else:
            truths = None
        index_dict = {}
        ### generate a dictionary that matches parameter names to column indices
        for index in range(4, len(header)):
            index_dict[header[index]] = index - 1
        lnL = samples[:,0]
        p = samples[:,1]
        p_s = samples[:,2]
        if cutoff != 0: # cutoff specified, so get the boolean mask
            mask = lnL > min_lnL
        elif frac != 1.0: # fraction specified but cutoff not, so get the appropriate mask
            ind = np.argsort(lnL)
            n = int(len(lnL) * frac)
            mask = ind[len(lnL) - n:]
        else: # no mask
            mask = [True] * len(lnL)
        lnL = lnL[mask]
        p = p[mask]
        p_s = p_s[mask]
        ### get columns of array corresponding to actual parameter samples
        x = samples[:,[index_dict[name] for name in params]]
        ### shift all the lnL values up so that we don't have rounding issues
        lnL -= np.max(lnL)
        L = np.exp(lnL)
        ### calculate weights
        weights = L * p / p_s
        print("eff samp:", np.sum(weights) / np.max(weights))
        weights /= np.sum(weights)
        ### throw out points with weight less than minimum weight
        mask2 = weights > min_weight
        print(np.sum(mask2), 'samples with weight >', min_weight)
        print('Median weight of full sample set:', np.median(weights))
        weights = weights[mask2]
        x = x[mask]
        x = x[mask2]
        print('Median weight of masked sample set:', np.median(weights))
        color = color_list[i % len(color_list)]
        for i in range(len(params)):
            if params[i] == 'log_mej_red':
                params[i] = 'mej_red'
                x[:,i] = 10.0**x[:,i]
            elif params[i] == 'log_mej_purple':
                params[i] = 'mej_purple'
                x[:,i] = 10.0**x[:,i]
            elif params[i] == 'log_mej_blue':
                params[i] = 'mej_blue'
                x[:,i] = 10.0**x[:,i]

        labels = []
        for param in params:
            if param in tex_dict:
                labels.append(tex_dict[param])
            else:
                labels.append(param)
        if combine and ind == len(total_samples) - 1:
            plot_density = True
        elif len(total_samples) == 1:
            plot_density = True
        else:
            plot_density = False
        ### make the corner plot
        fig_base = corner.corner(x, weights=weights, levels=args.cl, fig=fig_base, labels=labels, truths=truths,
                                 color=color, plot_datapoints=False, plot_density=plot_density,
                                 contours=True, smooth1d=0.1, smooth=0.1)
        i += 1
    if title is not None:
        plt.title(title)
    if leg is not None:
        ### figure out where to put the legend so it doesn't cover anything
        xcoord = 0 #len(params) - 1
        ycoord = len(params)
        ### generate the legend
        lgd = plt.legend(leg, bbox_to_anchor=(xcoord, ycoord), loc='upper left')
        #lgd = plt.legend(leg, loc="center left")
        ### fix the colors in the legend -- for some reason, if truth values are provided,
        ### every entry in the legend will have the same color
        for i in range(len(sample_files)):
            lgd.legendHandles[i].set_color(color_list[i])
        ### the extra arguments in savefig() make sure that the legend is not cut off
        plt.savefig(out, bbox_extra_artists=(lgd,), bbox_inches='tight')
    else:
        plt.savefig(out)

if __name__ == '__main__':
    args = _parse_command_line_args()
    sample_files = args.posterior_samples
    out = args.out
    params = args.p
    truths = args.truth_file
    cutoff = args.c
    frac = args.frac
    leg = args.legend
    cl = args.cl
    title = args.title
    combine = args.combine
    min_weight = args.min_weight
    generate_corner_plot(sample_files, out, params, truths, cutoff, frac, leg, cl,
                  title, combine, min_weight)
