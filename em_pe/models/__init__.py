from .interpolated_model import *
from .kilonova import *
from .kilonova_3c import *

model_dict = {
        "interpolated":interpolated,
        "kilonova":kilonova,
        "kilonova_3c":kilonova_3c
}

bounds_dict = {
    "dist":[10.0, 1000.0],
    "log_mej":[-3.0, -1.0],
    "vej":[0.1, 0.3],
    "log_kappa":[-1.0, 2.0],
    "Tc":[400.0, 5000.0],
    "log_mej_red":[-3.0, -1.0],
    "log_mej_purple":[-3.0, -1.0],
    "log_mej_blue":[-3.0, -1.0],
    "vej_red":[0.1, 0.3],
    "vej_purple":[0.1, 0.2],
    "vej_blue":[0.1, 0.2],
    "Tc_red":[400.0, 5000.0],
    "Tc_purple":[400.0, 5000.0],
    "Tc_blue":[400.0, 5000.0]
}
