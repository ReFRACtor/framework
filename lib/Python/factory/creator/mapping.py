from refractor import framework as rf

from .base import Creator
from .. import param


class Linear(Creator):
    '''Creator for Linear or "one-to-one" mapping'''

    def create(self, **kwargs):
        return rf.MappingLinear()


class Log(Creator):
    '''Creator for Log mapping'''

    def create(self, **kwargs):
        return rf.MappingLog()


class Offset(Creator):
    '''Creator for Offset mapping'''
    offset = param.Scalar(float, required=False, default=0.0)
    # TODO: Be nice to get offsetee from value() of our parent
    offsetee = param.Array(dims=1)

    def create(self, **kwargs):
        return rf.MappingOffset(self.offset(), self.offsetee())


class Scale(Creator):
    '''Creator for Scale mapping'''
    scaling = param.Scalar(float, required=False, default=1.0)
    # TODO: Be nice to get scalee from value() of our parent
    scalee = param.Array(dims=1)

    def create(self, **kwargs):
        return rf.MappingScale(self.scaling(), self.scalee())


class Gaussian(Creator):
    '''Creator for Gaussian mapping'''
    pressure = param.InstanceOf(rf.Pressure)
    linear_aod = param.Scalar(bool, default=True)

    def create(self, **kwargs):
        return rf.MappingGaussian(self.pressure(), self.linear_aod())

