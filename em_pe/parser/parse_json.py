# -*- coding: utf-8 -*-
'''
Parse JSON
----------
Parse Open Astronomy Catalog (OAC) JSON files and save requested data.

Example::

    $ python parse_json.py --f data.json --b R --b V --out Data/

To see full command line parameter documentation::

    $ python parse_json.py -h
    usage: parse_json.py [-h] [--t0 T0] [--f F] [--b B] [--out OUT]

    Parse Open Astronomy Catalog (OAC) JSON files

    optional arguments:
      -h, --help  show this help message and exit
      --t0 T0     Initial time (t=0 for event)
      --f F       Filename for JSON file
      --b B       Data bands to store
      --out OUT   Directory to save data to

'''

from __future__ import print_function
import numpy as np
import argparse
import sys
import json

def _parse_command_line_args():
    '''
    Parses and returns the command line arguments.
    '''
    parser = argparse.ArgumentParser(description='Parse Open Astronomy Catalog (OAC) JSON files')
    parser.add_argument('--t0', type=float, default=0, help='Initial time (t=0 for event)')
    parser.add_argument('--f', help='Filename for JSON file')
    parser.add_argument('--b', action='append', help='Data bands to store')
    parser.add_argument('--out', help='Directory to save data to')
    parser.add_argument('--maxpts', type=float, default=np.inf, help='Maximum number of points to keep for each band')
#    parser.add_argument('--tmax', type=float, default=np.inf, help='Upper bound for time points to keep')
#    parser.add_argument('--tmin', type=float, default=0, help='Lower bound for time points to keep')
    return parser.parse_args()

def _read_data(args):
    filename = args.f
    name = filename.split('/')[-1] # get rid of path except for filename
    name = name.split('.')[0] # get event name from filename
    ### read in the data
    with open(filename, "r") as read_file:
        data = json.load(read_file)[name]['photometry']
    ### create empty data arrays
    data_dict = {}
    for band in args.b:
        data_dict[band] = np.empty((4, 0))
    for entry in data:
        if 'band' in entry:
            band = entry['band']
            ### check that it's a band we want and that it has an error magnitude
            if band in args.b and 'e_magnitude' in entry:
                ### [time, time error, magnitude, magnitude error]
                to_append = np.array([[entry['time']], [0], [entry['magnitude']], [entry['e_magnitude']]]).astype(np.float)
                to_append[0] -= args.t0
                data_dict[band] = np.append(data_dict[band], to_append, axis=1)
    for band in data_dict:
        data = data_dict[band]
        ### check if we have too much data
        if data.shape[1] > args.maxpts:
            ### basically, generate random indices, take the columns (data points)
            ### specified by those columns, and then sort them based on times
            ### (sorting is not strictly necessary but it seems like a good idea
            ### to keep data ordered)
            cols = np.random.randint(0, data.shape[1], int(args.maxpts))
            data = data[:,cols]
            data = data[:,data[0].argsort()]
            data_dict[band] = data
    return data_dict

def _save_data(args, data_dict):
    outdir = args.out
    for band in data_dict:
        filename = outdir + band + '.txt'
        np.savetxt(filename, data_dict[band])

if __name__ == '__main__':
    args = _parse_command_line_args()
    data_dict = _read_data(args)
    _save_data(args, data_dict)
