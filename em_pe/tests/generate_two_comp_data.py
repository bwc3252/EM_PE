from __future__ import print_function
import numpy as np

from em_pe.models import two_comp
from em_pe.models.lightcurve_utils import calc_meje, calc_vej

np.random.seed(1)

dir = 'em_pe/tests/temp/' # directory to store data in

tini = 0.1
tmax = 8
n = 50 # number of samples to generate
tdays = np.linspace(tini, tmax, n)
err_lim = 0.05 # completely made-up based on typical errors seen in real data

mej1 = 0.01
vej1 = 0.1
mej2 = 0.03
vej2 = 0.2
frac = 0.5 # mixture ratio

print('mej1:', mej1, 'vej1:', vej1, 'mej2:', mej2, 'vej2:', vej2)

params = {'mej1':mej1, 'vej1':vej1, 'mej2':mej2, 'vej2':vej2, 'frac':frac}
t_bounds = [tini, tmax]

### initialize the model and set parameters

model = two_comp()
model.set_params(params, t_bounds)

### generate and save the data

for band in model.bands:
    filename = dir + band + '.txt'
    m, _ = model.evaluate(tdays, band)
    data = np.empty((4, n))
    data[0] = tdays
    data[2] = m + np.random.uniform(-1 * err_lim, err_lim, n) # generate errors
    data[3] = np.ones(n) * err_lim
    np.savetxt(filename, data.T)

### save true values

truths = np.array([mej1, vej1, mej2, vej2, frac])
filename = dir + 'test_truths.txt'
np.savetxt(filename, truths)
