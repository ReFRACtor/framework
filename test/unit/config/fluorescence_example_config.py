import numpy as np

from refractor.framework.factory import creator
from refractor.framework.config import refractor_config
from refractor import framework as rf

from .example_base_config import base_config

@refractor_config
def config(**kwargs):
    config_def = base_config()

    config_def['forward_model']['spectrum_effect']['fluorescence_effect'] = {
        'creator': creator.forward_model.FluorescenceEffect,
        'reference_point': {
            'creator': creator.value.ArrayWithUnit,
            'value': np.array([0.757]),
            'units':'micron',
        },
        'coefficients': np.array([-1.35039e-09, 0.0016]),
        'cov_unit': rf.Unit("ph / s / m^2 / micron sr^-1"),
        'which_channels': np.array([0]),
    }

    config_def['forward_model']['spectrum_effect']['effects'].append("fluorescence_effect")

    return config_def
