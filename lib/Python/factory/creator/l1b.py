import numpy as np

from .base import Creator
from .. import param

from refractor import framework as rf

class ValueFromLevel1b(Creator):
    "Extracts the surface pressure value from the attached meterological file"

    l1b = param.InstanceOf(rf.Level1b)
    field = param.Scalar(str)

    def create(self, channel_index=None, **kwargs):

        field_val = getattr(self.l1b(), self.field())
        if np.isscalar(field_val):
            return np.full(1, field_val)
        elif callable(field_val):
            if channel_index is not None:
                return field_val(channel_index)
            else:
                return field_val()
        else:
            return field_val
