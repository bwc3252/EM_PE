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
   low_level_utils

Examples and Usage
------------------

CLI Example for GW170817
^^^^^^^^^^^^^^^^^^^^^^^^
The following is an example of a bash script that will generate posterior
samples from GW170817 data (although these commands could be run directly from
the command line). A two-component model is used, and for simplicity, only two
bands are used::

    ### t=0 for event: 1187008882.43 for GW170817 (in GPS time)
    ### assumes temp/ directory exists
    python em_pe/parser/parse_json.py --t0 1187008882.43 --f Data/GW170817.json \
                            --b z --b y --out temp/ --maxpts 100 --tmax 8 \
                            --time_format gps

    n=60 # fix number of iterations to something reasonable

    ### assumes temp/z/ directory exists
    ### fix 'frac' parameter to 0.5 and 'dist' parameter to 40.0
    python em_pe/sampler.py --dat temp/ --m two_comp -v --f z.txt \
                            --min $n --max $n \
                            --out temp/z/posterior_samples.txt \
                            --fixed_param frac 0.5 --fixed_param dist 40 &

    ### assumes temp/y/ directory exists
    ### fix 'frac' parameter to 0.5 and 'dist' parameter to 40.0
    python em_pe/sampler.py --dat temp/ --m two_comp -v --f y.txt \
                            --min $n --max $n \
                            --out temp/y/posterior_samples.txt \
                            --fixed_param frac 0.5 --fixed_param dist 40

Once the posterior samples are generated, we can make a corner plot and a
lightcurve plot (again as a bash script)::

    python em_pe/plot_utils/plot_corner.py \
        --posterior_samples temp/z/posterior_samples.txt \
        --posterior_samples temp/y/posterior_samples.txt \
        --out temp/corner.png \
        --combine \
        --legend z \
        --legend y \
        --p mej1 --p vej1 --p mej2 --p vej2

    python -m em_pe/plot_utils/plot_lc --tmin 0.5 --tmax 8 --out temp/lc.png \
        --posterior_samples $dir/z/posterior_samples.txt \
        --posterior_samples $dir/y/posterior_samples.txt \
        --b z \
        --b y \
        --lc_file temp/z.txt \
        --lc_file temp/y.txt \
        --m two_comp --fixed_param frac 0.5 --fixed_param dist 40

Python API Example for GW170817
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The Python equivalent to the previous example::

    from em_pe import sampler
    from em_pe.parser import parse_json
    from em_pe.plot_utils import plot_corner, plot_lc

    t0 = 1187008882.43 # initial time for event
    data = 'Data/GW170817.json' # JSON file
    b = ['z', 'y'] # data bands to store
    dir = 'temp/' # directory to store data in

    ### Parse JSON data and convert to necessary format
    parse_json(t0, data, b, dir, maxpts=50, tmax=8, gps_time=True)

    model = 'two_comp'
    out = 'temp/z/posterior_samples.txt'
    n = 75
    fixed = [['frac', 0.5], ['dist', 40.0]]

    ### Initialize sampler (for z band)
    s = sampler(dir, model, ['z.txt'], out, min_iter=n, max_iter=n, fixed_params=fixed)

    ### Generate samples
    s.generate_samples()

    ### Initialize sampler (for y band)
    s = sampler(dir, model, ['y.txt'], out, min_iter=n, max_iter=n, fixed_params=fixed)

    ### Generate samples
    s.generate_samples()

    samples = ['temp/z/posterior_samples.txt', 'temp/y/posterior_samples.txt']
    out = temp/corner.png
    params = ['mej1', 'vej1', 'mej2', 'vej2']
    leg = ['z', 'y']

    ### Make a corner plot
    plot_corner.generate_corner_plot(samples, out, params, leg=leg)

    tmin = 0.5
    tmax = 8.0
    lc_file = ['temp/z.txt', 'temp/y.txt']

    ### Make a lightcurve plot
    plot_lc.generate_lc_plot(samples, out, model, tmin, tmax, b, lc_file=lc_file, fixed_params=fixed)
