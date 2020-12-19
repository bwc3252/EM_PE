# Sampler usage

Basic usage:

```bash
$ python3 em_pe/sampler.py --dat [directory] --m [model] --f [first data file] \
    --f [second data file] --min [min iterations] --max [max iterations] \
    --out [posterior sample filename]
```

Options:

- `--dat`: Directory containing the data files.
- `--m`: Name of model to use.
- `-v`: Verbose.
- `--cutoff`: Likelihood cutoff for storing posterior samples (default = 0).
- `--f`: Name of data file.
- `--min`: Minimum number of integrator iterations (default = 20).
- `--max`: Maximum number of integrator iterations (default = 20).
- `--ncomp`: Number of Gaussian components to use for a given parameter (e.g. `--ncomp mej 2`).
- `--fixed-param`: Fixed parameter name and value (e.g. `--fixed-param mej 0.01`).
- `--estimate-dist`: Option to fit distance.
- `--epoch`: Number of iterations before resetting sampling distributions.
- `--correlate-dims`: Parameters to group together for GMM sampler (e.g. `--correlate-dims mej vej`).
- `--burn-in`: Number of iterations for burn-in at start of sampling.
- `--beta-start`: Starting value for "beta" exponent used in burn-in.
- `--keep-npts`: Store the n highest-likelihood samples.
- `--nprocs`: Number of parallel processes to use for likelihood evaluations.
- `--set-limit`: Modify parameter limits (e.g. `--set limit mej 0.005 0.015`).
