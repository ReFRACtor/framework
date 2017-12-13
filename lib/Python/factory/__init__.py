from . import creator
from . import param

def process_config(config_def):

    if not isinstance(config_def, dict):
        raise TypeError("Configiguration definition is expected to be a dict instance")

    if "creator" in config_def:
        creator_class = config_def["creator"]
    else:
        creator_class = creator.base.SaveToCommon

    config_creator = creator_class(config_def)

    return config_creator()
