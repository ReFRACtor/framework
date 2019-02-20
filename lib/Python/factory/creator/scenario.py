import numpy as np

from .base import Creator
from .. import param
from .l1b import ValueFromLevel1b, RelativeAzimuthFromLevel1b

from refractor import framework as rf

# The 'scenario' values define the set of values that are needed by other creators
# that define the geometry, location and time of the radiances being modeled or sounding being retrieved
# 
# Below is the list of required values, each should be an ArrayWithUnit sized for the number of spectral channels of the instrument:
# time - list of Time objects 
# latitude - ArrayWithUnit(n_channels)
# longitude - ArrayWithUnit(n_channels)
# altitude - ArrayWithUnit(n_channels) 
# surface_height - ArrayWithUnit(n_channels) -- same as altitude, TBD to unify these two
# solar_zenith - ArrayWithUnit(n_channels)
# solar_azimuth - ArrayWithUnit(n_channels)
# observation_zenith - ArrayWithUnit(n_channels)
# observation_azimuth - ArrayWithUnit(n_channels)
# stokes_coefficient - array(n_channels, 4)

class ScenarioFromL1b(Creator):
    "Extracts necessary scenario values from a L1B file object"

    l1b = param.InstanceOf(rf.Level1b)

    def create(self, **kwargs):
        l1b_file = self.l1b()

        # Set up scenario values that have the same name as they are named in the L1B
        l1b_value_names = ['time', 'latitude', 'longitude', 'altitude', 'solar_zenith', 'solar_azimuth',
            "relative_velocity", "stokes_coefficient", "sample_grid"]
        values_from_l1b = { n:n for n in l1b_value_names }

        # Set up scenario values with different names
        values_from_l1b.update({
            "surface_height": "altitude",
            "observation_zenith": "sounding_zenith"
        })

        scenario_values = {}

        # Use ValueFromLevel1b logic to extract values
        for config_name, l1b_name in values_from_l1b.items():
            value_creator = ValueFromLevel1b({'l1b': l1b_file, 'field': l1b_name})
            scenario_values[config_name] = value_creator.create()

        scenario_values['observation_azimuth'] = RelativeAzimuthFromLevel1b({'l1b': l1b_file}).create()

        # Replicate values into the common store
        for config_name, value in scenario_values.items():
            self.common_store[config_name] = value

        return scenario_values
