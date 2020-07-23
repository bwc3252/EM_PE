.. EM_PE documentation master file, created by
   sphinx-quickstart on Tue Dec  4 23:24:43 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

EM PE Documentation
=================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   sampler
   plot_utils
   parser

Examples and Usage
------------------

Makefile Examples
^^^^^^^^^^^^^^^^^

A Makefile is provided to generate bash scripts for testing and common PE tasks.
For example, to generate fake data from the kilonova model and set up a PE run
to recover the parameters, run::

    $ make test_kilonova_3c

To run the parameter estimation and plotting codes::

    $ cd pe_runs/test_kilonova_3c*/
    $ ./sample.sh
    $ ./plot_corner.sh
    $ ./plot_lc.sh

The default settings assume `sample.sh` is being run on a machine with a large
number of CPU cores available, and as such will evaluate the likelihood in 8
parallel processes. The code runs on smaller machines as well, but the `--nprocs`
argument set in `sample.sh` should be set to something more appropriate.

The corner plot for this example should look similar to this:

[example image coming soon]

Using Log-Likelihood Function
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
In addition to full parameter estimation, this code also provides access to the
internal log-likelihood function::

    import numpy as np
    from em_pe import sampler

    dat = "./"
    m = "kilonova_3c"
    f = ["g.txt", ...]
    out = "placeholder.txt"

    ### Initialize sampler
    s = sampler(dat, m, f, out)

    ### Log-likelihood function takes a dictionary mapping parameter
    ### names to values

    params = {"mej_red":0.01, ...}

    ### Evaluate lnL

    lnL = s.log_likelihood(params)

    ### Alternatively, log_likelihood can take arrays of parameter samples

    n = 100

    params = {"mej_red":np.random.uniform(0.01, 0.04, n), ...}

    lnL = s.log_likelihood(params, vect=True)
