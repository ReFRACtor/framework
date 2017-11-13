import numpy as np

from .base import Creator
from .. import param

from refractor import framework as rf

class LidortRt(Creator):

    atmosphere = param.InstanceOf(rf.RtAtmosphere)
    stokes_coefficients = param.Array(dims=2)
    solar_zenith = param.ArrayWithUnit(dims=1)
    observation_zenith = param.ArrayWithUnit(dims=1)
    observation_azimuth = param.ArrayWithUnit(dims=1)

    num_streams = param.Scalar(int)
    num_mom = param.Scalar(int)

    pure_nadir = param.Scalar(bool, default=False)
    multiple_scattering_only = param.Scalar(bool, default=False)

    def create(self, **kwargs):
        stokes_object = rf.StokesCoefficientConstant(self.stokes_coefficients())

        return rf.LidortRt(self.atmosphere(), stokes_object, 
                self.solar_zenith().convert("deg").value, 
                self.observation_zenith().convert("deg").value, 
                self.observation_azimuth().convert("deg").value, 
                self.pure_nadir(), self.num_streams(), self.num_mom(), self.multiple_scattering_only())

class TwostreamRt(Creator):

    atmosphere = param.InstanceOf(rf.RtAtmosphere)
    stokes_coefficients = param.Array(dims=2)
    solar_zenith = param.ArrayWithUnit(dims=1)
    observation_zenith = param.ArrayWithUnit(dims=1)
    observation_azimuth = param.ArrayWithUnit(dims=1)

    do_fullquadrature = param.Scalar(bool, default=True)

    def create(self, **kwargs):
        stokes_object = rf.StokesCoefficientConstant(self.stokes_coefficients())

        
        return rf.TwostreamRt(self.atmosphere(), stokes_object,
                self.solar_zenith().convert("deg").value, 
                self.observation_zenith().convert("deg").value, 
                self.observation_azimuth().convert("deg").value, 
                self.do_fullquadrature())
