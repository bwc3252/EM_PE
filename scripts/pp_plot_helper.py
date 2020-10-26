import numpy as np
import argparse
import os


from em_pe.models import model_dict, param_dict

parser = argparse.ArgumentParser(description="Helper script to set up PP plot PE runs")
parser.add_argument("--m", help="Model to use")
parser.add_argument("--directory", help="Directory to set up PE runs in")
parser.add_argument("--name", help="Name to use for PP plot PE runs")
parser.add_argument("--npts", type=int, help="Number of different parameter sets to use")
parser.add_argument("--fixed-param", action="append", nargs=2, help="Set parameter with fixed value")
parser.add_argument("--sigma", default="0.0", help="Extra error estimate to be fit as a parameter")
parser.add_argument("--sampler-args", help="All extra arguments to pass to sampler in one string")
args = parser.parse_args()

base_dir = args.directory
if base_dir[-1] != "/":
    base_dir += "/"

base_dir += args.name + "/"
if not os.path.exists(base_dir):
    os.mkdir(base_dir)

m = model_dict[args.m]()
params = m.param_names + ["dist"]

fixed_params = {name:value for [name, value] in args.fixed_param} if args.fixed_param is not None else {}

base_args = "--n 100 --err 0.2 --time-format mjd --tmin 0.1 --tmax 30 --p sigma "  + args.sigma + " "
if args.sigma != "0.0":
    base_args += "--sigma " + args.sigma + " "
for p in fixed_params.keys():
    base_args += "--p " + p + " " + fixed_params[p] + " "

variable_params = {}
for p in params:
    if p not in fixed_params and p != "sigma":
        variable_params[p] = param_dict[p]()

commands = []

for i in range(args.npts):
    param_values = {}
    curr_dir = base_dir + str(i) + "/"
    if not os.path.exists(curr_dir):
        os.mkdir(curr_dir)
    command = "python3 ${EM_PE_INSTALL_DIR}/em_pe/tests/generate_data.py --m " + args.m + " --out " + str(i) + "/" + " " + base_args
    for p in variable_params.keys():
        param_values[p] = variable_params[p].sample_from_prior(width=1.0)
        command += " --p " + p + " " + str(param_values[p])
    commands.append(command)
    commands.append("python3 ${EM_PE_INSTALL_DIR}/em_pe/sampler.py --dat "
            + str(i) + "/" + " --m " + args.m + " -v --f g.txt --f r.txt --f "
            + "i.txt --f z.txt --f y.txt --f J.txt --f H.txt --f K.txt --min 40"
            + " --max 40 --out " + str(i) + "/" + "samples.txt"
            + " --burn-in 10 --beta-start 0.01 --keep-npts 100000 --nprocs 8 "
            + args.sampler_args)
    for p in fixed_params.keys():
        commands[-1] += " --fixed-param " + p + " " + str(fixed_params[p])
    for p in variable_params.keys():
        commands[-1] += " --set-limit " + p + " " + str(param_values[p] / 2.0) + " " + str(2.0 * param_values[p])

with open(base_dir + "run.sh", "w") as f:
    f.write("\n".join(commands))
