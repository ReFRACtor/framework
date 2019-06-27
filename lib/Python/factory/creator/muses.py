import numpy as np

from .base import Creator
from .. import param

from refractor import framework as rf
from . import atmosphere
from . import absorber

from refractor.muses_osp import OSP

class CreatorMUSES(Creator):

    l1b = param.InstanceOf(rf.Level1b)
    osp_directory = param.Scalar(str)

    def osp_common(self):

        l1b_obj = self.l1b()
        osp_dir = self.osp_directory()

        osp_common = {
            "base_dir": osp_dir,
            "latitude": l1b_obj.latitude(0).value,
            "longitude": l1b_obj.longitude(0).value,
            "obs_time": l1b_obj.time(0).as_datetime(),
        }

        return osp_common

    def pressure_obj(self, pressure_levels):
        return rf.PressureSigma(pressure_levels, pressure_levels[-1], False)

class PressureGridMUSES(CreatorMUSES, atmosphere.PressureGrid):

    num_levels = param.Scalar(int, default=48)
    pressure_levels = param.Choice(atmosphere.PressureGrid.pressure_levels, param.NoneValue())
    value = param.Choice(atmosphere.PressureGrid.value, param.NoneValue())

    def create(self, **kwargs):

        psur_osp = OSP("PSUR", use_log=False, **self.osp_common())

        if self.pressure_levels() is None:
            self.config_def["pressure_levels"] = psur_osp.pressure_resampled(self.num_levels())

        if self.value() is None:
            self.config_def["value"] = psur_osp.climatology_full_grid

        return atmosphere.PressureGrid.create(self, **kwargs)

class TemperatureMUSES(CreatorMUSES, atmosphere.TemperatureLevelOffset):

    temperature_levels = param.Choice(atmosphere.TemperatureLevelOffset.temperature_levels, param.NoneValue())
    pressure = param.Choice(atmosphere.TemperatureLevelOffset.pressure, param.NoneValue())

    def create(self, **kwargs):

        tatm_osp = OSP("TATM", use_log=False, **self.osp_common())
        
        if self.temperature_levels() is None:
            self.config_def["temperature_levels"] = tatm_osp.climatology_full_grid
            # Override pressure from that emitted through common_store
            self.config_def["pressure"] = self.pressure_obj(tatm_osp.pressure_full_grid)

        return atmosphere.TemperatureLevelOffset.create(self, **kwargs)

class SurfaceTemperatureMUSES(CreatorMUSES, atmosphere.SurfaceTemperature):

    num_channels = param.Scalar(int)
    value = param.Choice(atmosphere.SurfaceTemperature.value, param.NoneValue())

    def create(self, **kwargs):

        tsur_osp = OSP("TSUR", use_log=False, **self.osp_common())

        if self.value() is None:
            self.config_def['value'] = rf.ArrayWithUnit_double_1(np.full(self.num_channels(), tsur_osp.climatology_full_grid), "K")

        return atmosphere.SurfaceTemperature.create(self, **kwargs)

class AbsorberVmrMUSES(CreatorMUSES, absorber.AbsorberVmrLevel):

    def create(self, gas_name, **kwargs):
        gas_osp = OSP(gas_name, **self.osp_common())

        self.config_def['log_retrieval'] = True
        self.config_def['value'] = gas_osp.climatology_full_grid
        self.config_def['pressure'] = self.pressure_obj(gas_osp.pressure_full_grid)

        return absorber.AbsorberVmrLevel.create(self, gas_name, **kwargs)

class SurfacePressureCov(CreatorMUSES):

    def create(self, **kwargs):
        psur_osp = OSP("PSUR", use_log=False, **self.osp_common())
        return psur_osp.covariance_full_grid

class SurfaceTemperatureCov(CreatorMUSES):

    num_channels = param.Scalar(int)

    def create(self, **kwargs):
        tsur_osp = OSP("TSUR", use_log=False, **self.osp_common())
        return np.identity(self.num_channels()) * tsur_osp.covariance_full_grid

class AbsorberCov(CreatorMUSES):

    gas_name = param.Scalar(str)
    num_constraint_levels = param.Choice(param.Scalar(int), param.NoneValue(), required=False)

    def create(self, **kwargs):
        
        gas_osp = OSP(self.gas_name(), num_constraint_levels=self.num_constraint_levels(), **self.osp_common())

        # Not all gases have a covariance defined, those will not be retrieved
        if gas_osp.constraint_filename is not None:
            return np.linalg.inv(gas_osp.constraint_full_grid)

class AbsorberRetFlags(CreatorMUSES):

    gas_name = param.Scalar(str)
    num_constraint_levels = param.Choice(param.Scalar(int), param.NoneValue(), required=False)

    def create(self, **kwargs):
        
        gas_osp = OSP(self.gas_name(), num_constraint_levels=self.num_constraint_levels(), **self.osp_common())

        if gas_osp.constraint_filename is not None:
            return gas_osp.constraint_full_flags
        else:
            return np.ones(gas_osp.pressure_full_grid.shape[0], dtype=bool)

