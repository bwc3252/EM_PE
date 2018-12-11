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
        band = entry['band']
        ### check that it's a band we want and that it's not an "upper limit" 
        if band in args.b and ('upperlimit' not in entry or not entry['upperlimit']):
            ### [time, time error, magnitude, magnitude error]
            to_append = np.array([[entry['time']], [0], [entry['magnitude']], [entry['e_magnitude']]]).astype(np.float)
            to_append[0] -= args.t0
            data_dict[band] = np.append(data_dict[band], to_append, axis=1)
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
