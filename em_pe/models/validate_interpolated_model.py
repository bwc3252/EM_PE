# -*- coding: utf-8 -*-
"""
Validate Interpolated Models
----------------------------
"""
from __future__ import print_function
import argparse
import os
import numpy as np

parser = argparse.ArgumentParser(description="Check the accuracy of a previously interpolated model")
parser.add_argument("--m", default="", help="Name of model")
parser.add_argument("--tmin", type=float, default=0.0, help="Start time for interpolated model")
parser.add_argument("--tmax", type=float, default=10.0, help="Stop time for interpolated model")
parser.add_argument("--n", type=int, default=1000, help="Number of samples to use")
parser.add_argument("--n-times", type=int, default=50, help="Number of time points to use")
args = parser.parse_args()

if args.m == "me2017_lanthanide":
    from .me2017_lanthanide import me2017_lanthanide
