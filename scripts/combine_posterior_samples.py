# -*- coding: utf-8 -*-

import numpy as np
import argparse
import sys

parser = argparse.ArgumentParser(description='Combine multiple posterior sample files into a single file.')
parser.add_argument('--input-file', nargs="*", help='Input posterior sample file (can provide multiple instances)')
parser.add_argument('--output-fname', default='samples-combined.txt', help='Filename for output')
parser.add_argument('--keep-npts', type=int, help='Store the n highest-likelihood samples')
args = parser.parse_args()

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

if args.keep_npts is not None:
    ind_sorted = np.argsort(out[:,0])
    out = out[ind_sorted[out.shape[0] - args.keep_npts:]]

np.savetxt(args.output_fname, out, header=header)
