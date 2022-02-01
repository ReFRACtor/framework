from refractor.framework.factory import creator
from refractor.framework.config import refractor_config

from .example_base_config import base_config

@refractor_config
def config(**kwargs):
    config_def = base_config()

    # Use LRAD only
    config_def['radiative_transfer'] = {
        'creator': creator.rt.LRadRt,
        'num_streams': 4,
        'num_mom': 8,
    }

    # Use Rayleigh only.
    config_def['atmosphere']['aerosol'] = None

    return config_def
