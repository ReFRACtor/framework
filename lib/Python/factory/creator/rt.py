import numpy as np

from .base import Creator
from .. import param

from refractor import framework as rf

class LidortRt(Creator):

    atmosphere = param.InstanceOf(rf.RtAtmosphere)
    stokes_coefficient = param.Array(dims=2)
    solar_zenith = param.ArrayWithUnit(dims=1)
    observation_zenith = param.ArrayWithUnit(dims=1)
    observation_azimuth = param.ArrayWithUnit(dims=1)

    num_streams = param.Scalar(int)
    num_mom = param.Scalar(int)

    pure_nadir = param.Scalar(bool, default=False)
    multiple_scattering_only = param.Scalar(bool, default=False)

    use_solar_sources = param.Scalar(bool, default=True)
    use_thermal_emission = param.Scalar(bool, default=False)

    def create(self, **kwargs):
        stokes_object = rf.StokesCoefficientConstant(self.stokes_coefficient())

        return rf.LidortRt(self.atmosphere(), stokes_object, 
                self.solar_zenith().convert("deg").value, 
                self.observation_zenith().convert("deg").value, 
                self.observation_azimuth().convert("deg").value, 
                self.pure_nadir(), self.num_streams(), self.num_mom(), self.multiple_scattering_only(),
                self.use_solar_sources(), self.use_thermal_emission())

class TwostreamRt(Creator):

    atmosphere = param.InstanceOf(rf.RtAtmosphere)
    stokes_coefficient = param.Array(dims=2)
    solar_zenith = param.ArrayWithUnit(dims=1)
    observation_zenith = param.ArrayWithUnit(dims=1)
    observation_azimuth = param.ArrayWithUnit(dims=1)

    full_quadrature = param.Scalar(bool, default=True)

    use_solar_sources = param.Scalar(bool, default=True)
    use_thermal_emission = param.Scalar(bool, default=False)

    def create(self, **kwargs):
        stokes_object = rf.StokesCoefficientConstant(self.stokes_coefficient())
        
        return rf.TwostreamRt(self.atmosphere(), stokes_object,
                self.solar_zenith().convert("deg").value, 
                self.observation_zenith().convert("deg").value, 
                self.observation_azimuth().convert("deg").value, 
                self.full_quadrature(),
                self.use_solar_sources(), self.use_thermal_emission())

class LsiRt(Creator):

    atmosphere = param.InstanceOf(rf.RtAtmosphere)
    stokes_coefficient = param.Array(dims=2)
    solar_zenith = param.ArrayWithUnit(dims=1)
    observation_zenith = param.ArrayWithUnit(dims=1)
    observation_azimuth = param.ArrayWithUnit(dims=1)
    spec_win = param.InstanceOf(rf.SpectralWindow)

    lsi_config_file = param.Scalar(str)
    num_low_streams = param.Scalar(int)
    num_high_streams = param.Scalar(int)

    dedicated_twostream = param.Scalar(bool, default=True)
    pure_nadir = param.Scalar(bool, default=False)
    full_quadrature = param.Scalar(bool, default=True)
    use_lrad = param.Scalar(bool, default=True)
    
    def create(self, **kwargs):
        # Just use LIDORT for multiple scattering, when we use LRadRt
        do_multiple_scattering_only = self.use_lrad()

        stokes_object = rf.StokesCoefficientConstant(self.stokes_coefficient())

        # Minimum nmom allowed by LIDORT is 3
        nmom_low = min(self.num_low_streams() * 2, 3)

        if(self.num_low_streams() == 1 and self.dedicated_twostream()):
            rt_low = rf.TwostreamRt(self.atmosphere(), stokes_object,
                self.solar_zenith().convert("deg").value, 
                self.observation_zenith().convert("deg").value, 
                self.observation_azimuth().convert("deg").value, 
                self.full_quadrature())
        else:
            rt_low = rf.LidortRt(self.atmosphere(), stokes_object, 
                self.solar_zenith().convert("deg").value, 
                self.observation_zenith().convert("deg").value, 
                self.observation_azimuth().convert("deg").value, 
                self.pure_nadir(), self.num_low_streams(), nmom_low, do_multiple_scattering_only)

        lidort_pars = rf.Lidort_Pars.instance()
        rt_high = rf.LidortRt(self.atmosphere(), stokes_object, 
            self.solar_zenith().convert("deg").value, 
            self.observation_zenith().convert("deg").value, 
            self.observation_azimuth().convert("deg").value, 
            self.pure_nadir(), self.num_high_streams(), lidort_pars.maxmoments_input, do_multiple_scattering_only)

        spectral_bound = self.spec_win().spectral_bound

        if self.use_lrad():
            rt_low = rf.LRadRt(rt_low, spectral_bound, 
                self.solar_zenith().convert("deg").value, 
                self.observation_zenith().convert("deg").value, 
                self.observation_azimuth().convert("deg").value, 
                self.pure_nadir(), True, False)

            rt_high = rf.LRadRt(rt_high, spectral_bound, 
                self.solar_zenith().convert("deg").value, 
                self.observation_zenith().convert("deg").value, 
                self.observation_azimuth().convert("deg").value, 
                self.pure_nadir(), True, True)

        rt_high = rf.HresWrapper(rt_high)

        lsi_config = rf.HdfFile(self.lsi_config_file())
        return rf.LsiRt(rt_low, rt_high, lsi_config, "LSI")
