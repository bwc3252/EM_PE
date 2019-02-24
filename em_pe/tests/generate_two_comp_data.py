from __future__ import print_function
import numpy as np

from em_pe.models import two_comp
from em_pe.models.lightcurve_utils import calc_meje, calc_vej

np.random.seed(1)

tini = 0.1
tmax = 8
n = 50
tdays = np.linspace(tini, tmax, n)
err_lim = 0.05

mej1 = 0.01
vej1 = 0.1
mej2 = 0.03
vej2 = 0.2
frac = 0.5

print('mej1:', mej1, 'vej1:', vej1, 'mej2:', mej2, 'vej2:', vej2)

params = {'mej1':mej1, 'vej1':vej1, 'mej2':mej2, 'vej2':vej2, 'frac':frac}
t_bounds = [tini, tmax]

### initialize the model and set parameters

model = two_comp()
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

truths = np.array([mej1, vej1, mej2, vej2, frac])
filename = 'em_pe/tests/temp/test_truths.txt'
np.savetxt(filename, truths)
