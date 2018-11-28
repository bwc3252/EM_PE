from __future__ import print_function
import numpy as np
import argparse

from models import model_dict

def parse_command_line_args():
    parser = argparse.ArgumentParser(description='Generate posterior parameter samples from lightcurve data')
    parser.add_argument('--dat', help='Path to data directory')
    parser.add_argument('--m', action='append', help='Name of a model to use')
    parser.add_argument('--v', action='store_true', help='Verbose mode')
    parser.add_argument('--cutoff', default=0, type=float, help='Likelihood cutoff for storing posterior samples')
    parser.add_argument('--f', action='append', help='Name of a data file')
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_command_line_args()
    models = args.m
    data_files = [args.dat + file for file in args.f]
    print(models)
    print(data_files)
