from __future__ import print_function
import numpy as np

from em_pe.models import gwemlightcurve_models as models
from em_pe.models.lightcurve_utils import calc_meje, calc_vej

np.random.seed(1)

tini = 0.1
tmax = 10
n = 10
tdays = np.linspace(tini, tmax, n)
err_lim = 0.3
beta = 3.0
kappa_r = 1.0

m1 = 1.4
mb1 = 1.5
c1 = 0.15
m2 = 1.3
mb2 = 1.4
c2 = 0.2

mej = calc_meje(m1, mb1, c1, m2, mb2, c2)
vej = calc_vej(m1, c1, m2, c2)

params = {'mej':mej, 'vej':vej}
t_bounds = [tini, tmax]

### initialize the model and set parameters

model = models.woko2017()
model.set_params(params, t_bounds)

### generate and save the data

for band in model.bands:
    filename = 'em_pe/tests/temp/' + band + '.txt'
    m, _ = model.evaluate(tdays, band)
    data = np.empty((4, n))
    data[0] = tdays
    data[2] = m + np.random.uniform(-1 * err_lim, err_lim, n)
    data[3] = np.ones(n) * err_lim
    np.savetxt(filename, data)

### save true values

truths = np.array([mej, vej])
filename = 'em_pe/tests/temp/test_truths.txt'
np.savetxt(filename, truths)
