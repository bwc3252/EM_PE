# Code Integration Notes

Notes on integrating EM_PE code with [existing GW parameter estimation](https://github.com/oshaughn/research-projects-RIT/blob/temp-RIT-Tides-port_master-GPUIntegration/MonteCarloMarginalizeCode/Code/integrate_likelihood_extrinsic_batchmode).

## Overall structure of repository

```
EM_PE/

  Data/
    ##################################################
    # GW170817 data, model interpolation files, etc. #
    ##################################################

  doc/
    ########################################
    # documentation things for readthedocs #
    ########################################

  em_pe/
    integrator_utils/
      ########################################################
      # everything for the Monte Carlo integrator, GMM, etc. #
      ########################################################
      gaussian_mixture_model.py
      monte_carlo_integrator.py
      multivariate_truncnorm.py

    models/
      #####################################################################
      # everything for the lightcurve models                              #
      # NOTE: This directory's __init__.py file contains important things #
      #####################################################################
      lightcurve_utils.py
      linear_model.py           # simple test model
      me2017.py                 # currently has issues
      model.py                  # base class that all models inherit from
      test_models.py            # more simple test models
      two_comp_model.py         # two-component model based on woko2017
      woko2017.py

    parser/
      ############################################
      # currently just has the JSON parsing code #
      ############################################
      parse_json.py

    plot_utils/
      ##################################################
      # code to make lightcurve plots and corner plots #
      ##################################################
      plot_corner.py
      plot_lightcurves.py

    tests/
      ###################################################################
      # mostly just a place to dump temporary test data and put testing #
      # scripts                                                         #
      ###################################################################

    sampler.py # main user-facing code for generating posterior samples
```

## Usage

The basic workflow is as follows:

1. Parse JSON data:
```bash
$ python -m em_pe/parser/parse_json [options]
```
2. Generate posterior samples:
```bash
$ python -m em_pe/sampler [options]
```
3. Generate corner plot:
```bash
$ python -m em_pe/plot_utils/plot_corner [options]
```

Each of these steps can also be done using the Python API:

```Python
from em_pe import sampler
from em_pe.parser import parse_json
from em_pe.plot_utils import plot_corner

### Parse JSON data and save in correct format
parse_json.parse_json(...)

### Generate posterior samples
s = sampler.sampler(...)
s.generate_samples()

### Generate a corner plot
plot_corner.generate_plot(...)
```

## Options for Code Integration

I think there are three primary options for doing this:

1. Copy-paste this code into GW PE code
2. Implement or import GW likelihood in this sampling code
3. Import this code as a package and use Python API for EM likelihood

Options 2 and 3 make the most sense to me. If there is a straightforward way to
import and use the GW likelihood in this, that would be best:

```python
import gw_things

def log_likelihood(samples):
    lnL = ...
    lnL += gw_things.log_likelihood(samples)
    return lnL
```

This seems like the cleanest way to organize everything.

Otherwise, it would be possible to implement a Python API for the EM likelihood,
and something similar to the above code could be added to the GW likelihood.

## Questions

1. How hard would it be to import and use the GW likelihood?
2. Do we want to use my corner plot code, JSON parsing and associated internal
data format, and so on?
3. How important is it to support user-defined models? Currently, the user is
only able to choose between models that are already implemented. Are we assuming
that we will implement any models that someone might want to use?
