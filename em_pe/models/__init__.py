from .woko2017 import *
from .linear_model import *
from .woko2017_two_comp import *
from .oriented_model import *
from .interpolated_model import *
from .me2017_lanthanide import *
from .me2017_non_lanthanide import *
from .me2017_two_comp import *
from .test_model import *

model_dict = {
    'woko2017':woko2017,
    'linear':linear_model,
    'woko2017_two_comp':woko2017_two_comp,
    'oriented':oriented,
    'interpolated':interpolated,
    'me2017_lanthanide':me2017_lanthanide,
    'me2017_non_lanthanide':me2017_non_lanthanide,
    'me2017_two_comp':me2017_two_comp,
    'test':test
}

bounds_dict = {
    'a':[-1.0, 1.0],
    'b':[0.0, 1.0],
    'c':[-1.0, 1.0],
    'm':[0.0, 1.0],
    'y0':[20.0, 25.0],
    'mej':[1.0e-5, 0.5],
    'log_mej':[-16.0, -2.0],
    'mej_red':[1.0e-5, 0.5],
    'mej_blue':[1.0e-5, 0.5],
    'vej':[0.001, 0.9],
    'vej_red':[0.001, 0.9],
    'vej_blue':[0.001, 0.9],
    'frac':[0.0, 1.0],
    'dist':[10.0, 80.0],
    'angle':[-90.0, 90.0],
}
