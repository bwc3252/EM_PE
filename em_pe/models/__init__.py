from .test_models import *
from .me2017 import *
from .woko2017 import *
from .linear_model import *

model_dict = {
    'one_band_test':one_band_test_model,
    'two_band_test':two_band_test_model,
    'woko2017':woko2017,
    'me2017':me2017,
    'linear':linear_model
}

bounds_dict = {
    'a':[0.0, 2.0],
    'b':[0.0, 2.0],
    'm':[0.0, 1.0],
    'y0':[20.0, 25.0],
    'mej':[0.01, 0.1],
    'vej':[0.05, 0.2],
    'dist':[30.0, 50.0]
}
