from .woko2017 import *
from .woko2017_bns import *
from .interpolated_model import *
from .me2017_lanthanide import *
from .me2017_non_lanthanide import *

model_dict = {
        "woko2017":woko2017,
        "interpolated":interpolated,
        "me2017_lanthanide":me2017_lanthanide,
        "me2017_non_lanthanide":me2017_non_lanthanide,
        "woko2017_bns":woko2017_bns
}

bounds_dict = {
    'mej':[1.0e-5, 0.3],
    'vej':[0.001, 0.5],
    'dist':[10.0, 1000.0],
    'm1':[1.1, 1.90],
    'm2':[1.1, 1.90]
}
