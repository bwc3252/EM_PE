from .woko2017 import *
from .linear_model import *
from .two_comp_model import *
from .oriented_model import *
from .interpolated_model import *
from .me2017_lanthanide import *
from .me2017_non_lanthanide import *

model_dict = {
    'woko2017':woko2017,
    'linear':linear_model,
    'two_comp':two_comp,
    'oriented':oriented,
    'interpolated':interpolated,
    'me2017_lanthanide':me2017_lanthanide,
    'me2017_non_lanthanide':me2017_non_lanthanide
}

bounds_dict = {
    'a':[-1.0, 1.0],
    'b':[-1.0, 1.0],
    'c':[-1.0, 1.0],
    'm':[0.0, 1.0],
    'y0':[20.0, 25.0],
    'mej':[1.0e-7, 1.0],
    'vej':[0.01, 0.99],
    'mej1':[1.0e-7, 1.0],
    'mej2':[1.0e-7, 1.0],
    'vej1':[0.01, 0.99],
    'vej2':[0.01, 0.99],
    'frac':[0.0, 1.0],
    'dist':[30.0, 50.0],
    'angle':[0.0, 90.0],
    'kappa_r_nl':[0.7, 5.0], # non-lanthanide
    'kappa_r_l':[5.0, 15.0] # lanthanide
}
