import numpy as np

import refractor.factory.creator as creator
from refractor.config import refractor_config
from refractor import framework as rf

from .example_base_config import base_config, num_channels

@refractor_config
def spectral_fit_config(**kwargs):

    config_def = base_config(**kwargs)

    config_def['instrument']['dispersion'] = {
        'creator': creator.instrument.SpectralFitSampleGridCreator,
        'spectral_domains': {
            'creator': creator.l1b.ValueFromLevel1b,
            'field': 'sample_grid',
        },
        # Size of shift, squeeze the number of sensors (channels)
        # Create an array with the same value for all sensors
        'shift': np.full(num_channels, 1),
        'squeeze': np.full(num_channels, 0),
    }

    return config_def
