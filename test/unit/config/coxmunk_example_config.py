from refractor.framework.config import refractor_config

from .example_base_config import base_config

@refractor_config
def config(**kwargs):
    config_def = base_config()

    config_def['atmosphere']['ground']['child'] = 'coxmunk'

    return config_def
