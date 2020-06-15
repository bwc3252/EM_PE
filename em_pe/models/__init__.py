from .interpolated_model import *
from .kilonova import *
from .kilonova_2c import *
from .kilonova_3c import *

model_dict = {
        "interpolated":interpolated,
        "kilonova":kilonova,
        "kilonova_2c":kilonova_2c,
        "kilonova_3c":kilonova_3c
}

bounds_dict = {
    "dist":[35.0, 45.0],
    "log_mej":[-3.0, -1.0],
    "vej":[0.1, 0.4],
    "log_kappa":[-1.0, 2.0],
    "Tc":[400.0, 5000.0],
    "log_mej_red":[-3.0, -1.0],
    "log_mej_purple":[-3.0, -1.0],
    "log_mej_blue":[-3.0, -1.0],
    "vej_red":[0.1, 0.4],
    "vej_purple":[0.1, 0.4],
    "vej_blue":[0.1, 0.4],
    "Tc_red":[400.0, 5000.0],
    "Tc_purple":[400.0, 5000.0],
    "Tc_blue":[400.0, 5000.0],
    "t_exp":[-0.5, 0.5]
}
