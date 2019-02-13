from __future__ import print_function
import numpy as np

from em_pe.models import woko2017
from em_pe.models.lightcurve_utils import calc_meje, calc_vej

np.random.seed(1)

tini = 0.1
tmax = 8
n = 50
tdays = np.linspace(tini, tmax, n)
err_lim = 0.05

mej = 0.01 #calc_meje(m1, mb1, c1, m2, mb2, c2)
vej = 0.1 #calc_vej(m1, c1, m2, c2)

print('mej:', mej, 'vej:', vej)

params = {'mej':mej, 'vej':vej}
t_bounds = [tini, tmax]

### initialize the model and set parameters

model = woko2017()
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
