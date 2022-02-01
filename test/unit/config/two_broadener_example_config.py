from refractor.framework.factory import creator
from refractor.framework.config import refractor_config

from .example_base_config import base_config, static_value

@refractor_config
def config(**kwargs):
    config_def = base_config()

    config_def['atmosphere']['absorber']['gases'].append("CF4")

    # This ABSCO file has 2 broadeners
    config_def['atmosphere']['absorber']['O2']['absorption']['creator'] = creator.absorber.AbscoAer
    config_def['atmosphere']['absorber']['O2']['absorption']['filename'] = "{absco_base_path}/aer/O2_06200-06201_v0.0_init.nc"
    config_def['atmosphere']['absorber']['O2']['absorption']['interp_method'] = creator.AbscoInterpolationOption.nearest_neighbor

    # This ABSCO file has no broadeners
    config_def['atmosphere']['absorber']['CF4'] = {
        'creator': creator.absorber.AbsorberGasDefinition,
        'vmr': {
            'creator': creator.absorber.AbsorberVmrLevel,
            'vmr_profile': {
                'creator': creator.atmosphere.ConstantForAllLevels,
                # Nonsense value, we just need something so we copy the O2 apriori
                'value': static_value("Gas/O2/average_mole_fraction")[0],
            },
            'mapping': creator.mapping.NotRetrieved,
        },
        'absorption': {
            'creator': creator.absorber.AbscoAer,
            'table_scale': 1.0,
            'interp_method': creator.AbscoInterpolationOption.nearest_neighbor,
            'filename': "{absco_base_path}/aer/CH4_06200-06201_v0.0_init.nc",
        },
    }

    return config_def
