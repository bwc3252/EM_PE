from .test_models import *
from .me2017 import *
from .woko2017 import *
from .linear_model import *
from .two_comp_model import *
from .oriented_model import *
from .lightcurve_utils import *

model_dict = {
    'one_band_test':one_band_test_model,
    'two_band_test':two_band_test_model,
    'woko2017':woko2017,
    'me2017':me2017,
    'linear':linear_model,
    'two_comp':two_comp,
    'oriented':oriented
}

bounds_dict = {
    'a':[0.0, 2.0],
    'b':[0.0, 2.0],
    'm':[0.0, 1.0],
    'y0':[20.0, 25.0],
    'mej':[0.001, 0.1],
    'vej':[0.01, 0.5],
    'mej1':[0.001, 0.1],
    'vej1':[0.01, 0.5],
    'mej2':[0.001, 0.1],
    'vej2':[0.01, 0.5],
    'frac':[0.0, 1.0],
    'dist':[30.0, 50.0],
    'angle':[0.0, 90.0]
}
