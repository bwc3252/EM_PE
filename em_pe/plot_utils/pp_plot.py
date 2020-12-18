import numpy as np
from scipy.stats import norm
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import argparse
import os
import json

from em_pe.models import model_dict, param_dict

parser = argparse.ArgumentParser(description="Script to generate PP plot")
parser.add_argument("--m", help="Model to use")
parser.add_argument("--directory", help="PP plot run directory")
parser.add_argument("--name", help="Name")
parser.add_argument("--exclude-param", action="append", help="parameter to exclude from plot")
parser.add_argument("--exclude-dir", action="append", help="Exclude a run by directory name")
args = parser.parse_args()

exclude_params = args.exclude_param if args.exclude_param is not None else []
exclude_dir = args.exclude_dir if args.exclude_dir is not None else []

eff_samp_cutoff = 10.0

def binomial_credible_interval(phat, n, z):
    offset = np.sqrt(phat * (1.0 - phat) / n + z**2 / (4 * n**2))
    my_1 = np.ones((len(phat), 2))
    return 1.0 / (1.0 + z**2 / n) *  (np.outer(phat, np.array([1, 1])) + my_1 * z**2 / (2 * n) + np.outer(z * offset, np.array([-1, 1])))
    
pvalFiducial = 0.9 # fiducial credible level for outermost CI.  (Will change depending on number of parameters)

def binomial_credible_interval_default(phat, n, nParams=2):
    z_fiducial = norm.ppf((np.power(1.0 - (1.0 - pvalFiducial) / 2.0 , 1.0 / nParams)))
    #print("fiducial z", z_fiducial, " for nParams = ", nParams)
    return binomial_credible_interval(phat, n, z_fiducial)

tex_dict = {'mej':'$m_{ej}$',
            'mej_red':'$m_{ej}$ [red]',
            'mej_purple':'$m_{ej}$ [purple]',
            'mej_blue':'$m_{ej}$ [blue]',
            'vej':'$v_{ej}$',
            'vej_red':'$v_{ej}$ [red]',
            'vej_purple':'$v_{ej}$ [purple]',
            'vej_blue':'$v_{ej}$ [blue]',
            'Tc_red':'Tc [red]',
            'Tc_purple':'Tc [purple]',
            'Tc_blue':'Tc [blue]',
            'dist':'Distance',
            'kappa':'$\kappa$',
            'sigma':'$\\sigma$'
}

base_dir = args.directory
if base_dir[-1] != "/":
    base_dir += "/"

base_dir += args.name
if base_dir[-1] != "/":
    base_dir += "/"

m = model_dict[args.m]()
ordered_params = m.param_names + ["dist"]
params = {p:param_dict[p]() for p in ordered_params}
params_used = set()

cache = {}
eff_samp_list = []
dir_list = []

for d in os.listdir(base_dir):
    curr_dir = base_dir + d + "/"
    if not os.path.isdir(curr_dir) or "samples.txt" not in os.listdir(curr_dir) or curr_dir in exclude_dir or d in exclude_dir:
        continue
    print("reading samples from", curr_dir)
    s = np.loadtxt(curr_dir + "samples.txt")
    truths = np.loadtxt(curr_dir + "test_truths.txt")
    truths = dict(zip(ordered_params, truths))
    with open(curr_dir + "samples.txt") as f:
        ### the "header" contains the column names
        header = f.readline().strip().split(" ")[1:]
    samples = {header[i]:s[:,i] for i in range(len(header))}
    lnL = samples["lnL"]
    p = samples["p"]
    p_s = samples["p_s"]
    temp = np.exp(lnL) * p / p_s
    eff_samp = np.sum(temp) / np.max(temp)
    print("    eff_samp = {}".format(eff_samp))
    eff_samp_list.append(eff_samp)
    dir_list.append(d)
    L = np.exp(lnL - np.max(lnL))
    weights = L * p / p_s
    weights /= np.sum(weights)
    for i in range(3, len(header)):
        p = header[i]
        if p in exclude_params: continue
        if p not in params_used:
            cache[p] = {"CDF":[], "true":[], "ML":[]}
            params_used.add(p)
        cdf = np.sum(weights[samples[p] < truths[p]])
        cache[p]["CDF"].append(cdf)
        cache[p]["true"].append(truths[p])
        cache[p]["ML"].append(samples[p][np.argmax(lnL)])

eff_samp = np.array(eff_samp_list)
mask = eff_samp > eff_samp_cutoff
npts = eff_samp.size
npts_valid = np.sum(mask)

out_dict = {
        i:{"params":{p:{} for p in params_used}, "eff_samp":eff_samp[i], "directory":dir_list[i]
    } for i in range(npts)
}

plt.figure(figsize=(8, 8))

markers = {
        "mej":"x",
        "mej_red":"x",
        "mej_purple":"x",
        "mej_blue":"x",
        "vej":"o",
        "vej_red":"o",
        "vej_purple":"o",
        "vej_blue":"o",
        "sigma":"+"
}
colors = {
        "mej":"red",
        "mej_red":"red",
        "mej_purple":"purple",
        "mej_blue":"blue",
        "vej":"blue",
        "vej_red":"red",
        "vej_purple":"purple",
        "vej_blue":"blue",
        "sigma":"black"
}

for p in params_used:
    if p in exclude_params: continue
    cdf = np.array(cache[p]["CDF"])
    ind = np.argsort(cdf[mask])
    j = 0.0
    for i in range(npts):
        out_dict[i]["params"][p]["CDF"] = cache[p]["CDF"][i]
        out_dict[i]["params"][p]["true"] = cache[p]["true"][i]
        out_dict[i]["params"][p]["ML"] = cache[p]["ML"][i]
        if mask[i]:
            (xval,) = np.where(ind == j)[0] / npts_valid
            j += 1.0
        else:
            xval = np.nan
        out_dict[i]["params"][p]["xval"] = xval
    cdf = cdf[mask][ind]
    plt.scatter(np.arange(npts_valid) / float(npts_valid), cdf,
            label=(tex_dict[p] if p in tex_dict.keys() else p),
            marker=(markers[p] if p in markers.keys() else None),
            color=(colors[p] if p in colors.keys() else None)
    )
    


xvals = np.linspace(0.0, 1.0, 100)
plt.plot(xvals, xvals, "--", color="black")
pvals_lims = binomial_credible_interval_default(xvals, npts_valid, nParams=(len(params_used) + len(exclude_params)))
plt.plot(xvals, pvals_lims[:,0], color='k', ls=':')
plt.plot(xvals, pvals_lims[:,1], color='k', ls=':')

plt.xlabel("$P(x_{\\rm inj})$")
plt.ylabel("$\hat{P}$")
plt.legend()
plt.tight_layout()

plt.savefig("pp_plot.png")

with open("out.json", "w") as f:
    json.dump(out_dict, f, indent=4)
