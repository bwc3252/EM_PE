from __future__ import print_function
import numpy as np

a = 1
b = 1
err_lim = 0.01

tmin = -4
tmax = 4
n = 50
t = np.linspace(tmin, tmax, n)

data = np.empty((4, n))
data[0] = t
data[2] = (1 / np.cosh(a * t)) + b + np.random.uniform(-1 * err_lim, err_lim, n)
data[3] = np.ones(n) * err_lim

filename = 'em_pe/tests/Data/test_bandA.txt'
np.savetxt(filename, data)

truths = np.array([a, b])
filename = 'em_pe/tests/Data/test_truths.txt'
np.savetxt(filename, truths)
