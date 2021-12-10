from .base import Creator
from .. import param

from .atmosphere import AtmosphereCreator
from .forward_model import ForwardModel

import refractor.framework as rf

class PressureWithCloudHandling(Creator):
    """Creates a PressureWithCloudHandling object that encapsulates another Pressure object 
    with the ability to switch on and off cloud computations when used in conjunction with
    ForwardModelWithCloudHandling"""

    pressure_clear = param.InstanceOf(rf.Pressure)
    cloud_pressure_level = param.Scalar(float)

    def create(self, **kwargs):

        return rf.PressureWithCloudHandling(self.pressure_clear(), self.cloud_pressure_level())

class GroundWithCloudHandling(Creator):
    """Creates a GroundWithCloudHandling object that encapsulates another Ground object 
    with the ability to switch on and off cloud computations when used in conjunction with
    ForwardModelWithCloudHandling"""

    ground_clear = param.InstanceOf(rf.Ground)
    cloud_albedo = param.Scalar(float)

    def create(self, **kwargs):

        return rf.GroundWithCloudHandling(self.ground_clear(), self.cloud_albedo())

class AtmosphereWithCloudHandling(AtmosphereCreator):

    """This creator stands in place for the nominal atmosphere creator. It wraps the Pressure and Ground objects 
    with those that enable the ability to model clouds using these values:

    * cloud_pressure_level
    * cloud_albedo

    This creator is convenience over using PressureWithCloudHandling and GroundWithCloudHandling seperately
    by replacing the pressure and ground values with in config definition

    This creator must be used with the ForwardModelWithCloudHandling creator which defines:
    * cloud_fraction
    """

    cloud_pressure_level = param.Scalar(float)
    cloud_albedo = param.Scalar(float)

    def __init__(self, config_def, common_store=None):

        # Proxy the 
        config_def['pressure'] = { 
            'creator': PressureWithCloudHandling,
            'pressure_clear': config_def['pressure'],
            'cloud_pressure_level': config_def['cloud_pressure_level'],
        }

        config_def['ground'] = {
            'creator': GroundWithCloudHandling,
            'ground_clear': config_def['ground'],
            'cloud_albedo': config_def['cloud_albedo'],
        }
 
        super().__init__(config_def, common_store)

class ForwardModelWithCloudHandling(ForwardModel):

    """Creates a ForwardModelWithCloudHandling object that handles computing cloudy atmospheres defined by:
    * cloud_fraction

    It must be used with a AtmosphereWithCloudHandling creator or with the PressureWithCloudHandling and
    GroundWithCloudHandling creators together which define:
    * cloud_pressure_level
    * cloud_albedo

    Currently, the only retrievable component is cloud_fraction.
    """

    cloud_fraction = param.Scalar(float)

    def create(self, **kwargs):

        fm_clear = super().create(**kwargs)
        cloud_fraction = rf.CloudFractionFromState(self.cloud_fraction())

        fm_cloud = rf.ForwardModelWithCloudHandling(fm_clear, cloud_fraction)

        # Monkey patch this attribute since it is needed in other Creators
        # TODO: Add spectral_grid to forward model interface or move out to its own Creator
        fm_cloud.spectral_grid = fm_clear.spectral_grid

        return fm_cloud
