from refractor import framework as rf

from .base import Creator
from .. import param


class Linear(Creator):
    '''Creator for Linear or "one-to-one" mapping'''

    def create(self, **kwargs):
        return rf.StateMappingLinear()


class Log(Creator):
    '''Creator for Log mapping'''

    def create(self, **kwargs):
        return rf.StateMappingLog()


class Offset(Creator):
    '''Creator for Offset mapping'''
    offset = param.Scalar(float, required=False, default=0.0)
    # TODO: Be nice to get offsetee from value() of our parent
    offsetee = param.Array(dims=1)

    def create(self, **kwargs):
        return rf.StateMappingOffset(self.offset(), self.offsetee())


class Scale(Creator):
    '''Creator for Scale mapping'''
    scaling = param.Scalar(float, required=False, default=1.0)
    # TODO: Be nice to get scalee from value() of our parent
    scalee = param.Array(dims=1)

    def create(self, **kwargs):
        return rf.StateMappingScale(self.scaling(), self.scalee())


class Gaussian(Creator):
    '''Creator for Gaussian mapping'''
    pressure = param.InstanceOf(rf.Pressure)
    linear_aod = param.Scalar(bool, default=True)

    def create(self, **kwargs):
        return rf.StateMappingGaussian(self.pressure(), self.linear_aod())


class InterpolateLinearLinear(Creator):
    '''Creator for Gaussian mapping'''
    pressure_to = param.InstanceOf(rf.Pressure)
    pressure_from = param.InstanceOf(rf.Pressure)

    def create(self, **kwargs):
        return rf.StateMappingInerpolateLinearLinear(self.pressure_to(), self.pressure_from())

class InterpolateLogLinear(Creator):
    '''Creator for Gaussian mapping'''
    pressure_to = param.InstanceOf(rf.Pressure)
    pressure_from = param.InstanceOf(rf.Pressure)

    def create(self, **kwargs):
        return rf.StateMappingInerpolateLogLinear(self.pressure_to(), self.pressure_from())

class InterpolateLogLog(Creator):
    '''Creator for Gaussian mapping'''
    pressure_to = param.InstanceOf(rf.Pressure)
    pressure_from = param.InstanceOf(rf.Pressure)

    def create(self, **kwargs):
        return rf.StateMappingInerpolateLogLog(self.pressure_to(), self.pressure_from())

class InterpolateLinearLog(Creator):
    '''Creator for Gaussian mapping'''
    pressure_to = param.InstanceOf(rf.Pressure)
    pressure_from = param.InstanceOf(rf.Pressure)

    def create(self, **kwargs):
        return rf.StateMappingInerpolateLinearLog(self.pressure_to(), self.pressure_from())

class Composite(Creator):
    '''Creator for Gaussian mapping'''

    mappings = param.Iterable(rf.StateMapping)

    def create(self, **kwargs):

        mappings_inp = self.mappings()
        mappings_vec = rf.vector_state_mapping()

        for map_obj in mappings_inp:
            mappings_vec.push_back(map_obj)
        
        return rf.StateMappingComposite(mappings_vec)
