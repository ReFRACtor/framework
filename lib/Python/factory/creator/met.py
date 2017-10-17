import numpy as np

from .base import Creator
from .. import param

from refractor import framework as rf

class ValueFromMet(Creator):
    "Extracts the surface pressure value from the attached meterological file"

    met = param.InstanceOf(rf.Meteorology)
    field = param.Scalar(str)

    def create(self, **kwargs):

        field_val = getattr(self.met(), self.field())
        if np.isscalar(field_val):
            return np.full(1, field_val)
        else:
            return field_val
