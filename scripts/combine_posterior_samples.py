# -*- coding: utf-8 -*-

import numpy as np
import argparse
import sys

parser = argparse.ArgumentParser(description='Combine multiple posterior sample files into a single file.')
parser.add_argument('--input-file', nargs="*", help='Input posterior sample file (can provide multiple instances)')
parser.add_argument('--output-fname', default='samples-combined.txt', help='Filename for output')
parser.add_argument('--keep-npts', type=int, help='Store the n highest-likelihood samples')
parser.add_argument('--tempering-exp', default=1.0, type=float, help="Exponent for likelihoods")
parser.add_argument('--max-lnL', default=np.inf, type=float, help="Maximum log-likelihood")
parser.add_argument('--set-limit', action='append', nargs=3, help='Modify parameter limits (e.g. --set-limit mej_red 0.008 0.012)')
args = parser.parse_args()
    
if args.set_limit is not None:
    limits = {name:(float(llim), float(rlim)) for (name, llim, rlim) in args.set_limit}
else:
    limits = None

if args.input_file is None:
    print("No input files, exiting")
    sys.exit()

out = []
for fname in args.input_file:
    print("Loading samples from {}...".format(fname))
    out.append(np.loadtxt(fname))

out = np.concatenate(out, axis=0)

with open(args.input_file[0], "r") as f:
    header = f.readline()[2:] # read the header, remove the "# " at the beginning
    header = header[:-1] # remove the '\n' at the end

out = out[out[:,0] < args.max_lnL]
out[:,0] *= args.tempering_exp

if args.keep_npts is not None:
    ind_sorted = np.argsort(out[:,0])
    out = out[ind_sorted[out.shape[0] - args.keep_npts:]]

if limits is not None:
    params = header.split(" ")
    for i, p in enumerate(params):
        if p in limits.keys():
            llim, rlim = limits[p]
            out = out[out[:,i] > llim]
            out = out[out[:,i] < rlim]


np.savetxt(args.output_fname, out, header=header)
