import numpy as np
from scipy.interpolate import interp1d
import argparse

parser = argparse.ArgumentParser(description="Script to generate injection/recovery tests for lightcurves extracted from simulations")
parser.add_argument("--input-sim", help="File containing simulation magnitudes")
parser.add_argument("--tmin", type=float, help="Minimum time value")
parser.add_argument("--tmax", type=float, help="Maximum time value")
parser.add_argument("--n", type=int, help="Number of points per band")
parser.add_argument("--err", type=float, help="Error std. dev.")
parser.add_argument("--out", help="Directory to write output to")
parser.add_argument("--angular-bin", type=int, help="Angular bin from which to pull data")
parser.add_argument("--dist", type=float, help="Distance in Mpc")
args = parser.parse_args()

bands_out = ["g", "r", "i", "z", "y", "J", "H", "K"]

# this line loads the data, splits it by band, removes the last one (the S band), and makes a dictionary mapping bands to data
data_in = dict(zip(bands_out, np.split(np.loadtxt(args.input_sim), 9)[:-1]))

# grab the time values
t_in = data_in["g"][:,1]

# remove all the data except the specific mags we want
data_in = {b:data_in[b][:,args.angular_bin + 2] for b in bands_out}

for b in bands_out:
    f = interp1d(t_in, data_in[b], fill_value="extrapolate")
    t = np.exp(np.random.uniform(np.log(args.tmin), np.log(args.tmax), args.n))
    t = np.sort(t)
    out = np.empty((args.n, 4))
    out[:,0] = t
    out[:,2] = f(t) + np.random.normal(0.0, args.err, args.n) + 5.0 * (np.log10(args.dist * 1.0e6) - 1.0)
    out[:,3] = args.err
    np.savetxt(args.out + b + ".txt", out)
