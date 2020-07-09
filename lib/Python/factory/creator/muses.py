import re
import numpy as np

from .base import Creator
from .. import param

from refractor import framework as rf
from . import atmosphere
from . import absorber
from . import retrieval

from refractor.muses_osp import OSP

class CreatorMUSES(Creator):

    l1b = param.InstanceOf(rf.Level1b)
    osp_directory = param.Scalar(str)
    osp_instrument_name = param.Scalar(str, default=None)

    def osp_common(self):

        l1b_obj = self.l1b()

        osp_common = {
            "base_dir": self.osp_directory(),
            "latitude": l1b_obj.latitude(0).value,
            "longitude": l1b_obj.longitude(0).value,
            "obs_time": l1b_obj.time(0).as_datetime(),
            "instrument_name": self.osp_instrument_name(),
        }

        return osp_common

    def pressure_obj(self, pressure_levels):
        return rf.PressureSigma(pressure_levels, pressure_levels[-1], False)

class PressureGridMUSES(CreatorMUSES, atmosphere.PressureGrid):

    num_levels = param.Scalar(int, default=48)
    pressure_levels = param.Choice(atmosphere.PressureGrid.pressure_levels, param.NoneValue())
    value = param.Choice(atmosphere.PressureGrid.value, param.NoneValue())

    def create(self, **kwargs):

        psur_osp = OSP("PSUR", **self.osp_common())

        if self.pressure_levels() is None:
            self.config_def["pressure_levels"] = psur_osp.pressure_resampled(self.num_levels())

        if self.value() is None:
            self.config_def["value"] = psur_osp.climatology_full_grid

        return atmosphere.PressureGrid.create(self, **kwargs)

class TemperatureMUSES(CreatorMUSES, atmosphere.TemperatureLevel):

    value = param.Choice(atmosphere.TemperatureLevel.value, param.NoneValue())

    def create(self, **kwargs):

        tatm_osp = OSP("TATM", **self.osp_common())

        ret_levels = tatm_osp.retrieval_levels
        
        if self.value() is None:
            self.config_def["value"] = tatm_osp.climatology_full_grid[ret_levels]
        else:
            inp_temperature = self.value()
            inp_pressure = self.pressure().pressure_grid.value.value
            resamp_temp = np.interp(tatm_osp.pressure_full_grid, inp_pressure, inp_temperature)

            self.config_def["value"] = resamp_temp[ret_levels]

        self.config_def["pressure"] = self.pressure_obj(tatm_osp.pressure_full_grid[ret_levels])

        return atmosphere.TemperatureLevel.create(self, **kwargs)

class SurfaceTemperatureMUSES(CreatorMUSES, atmosphere.SurfaceTemperature):

    num_channels = param.Scalar(int)
    value = param.Choice(atmosphere.SurfaceTemperature.value, param.NoneValue())

    def create(self, **kwargs):

        tsur_osp = OSP("TSUR", **self.osp_common())

        if self.value() is None:
            self.config_def['value'] = rf.ArrayWithUnit_double_1(np.full(self.num_channels(), tsur_osp.climatology_full_grid), "K")

        return atmosphere.SurfaceTemperature.create(self, **kwargs)

class AbsorberVmrMUSES(CreatorMUSES, absorber.AbsorberVmrLevel):

    def create(self, gas_name, **kwargs):
        gas_osp = OSP(gas_name, **self.osp_common())

        ret_levels = gas_osp.retrieval_levels

        self.config_def['log_retrieval'] = gas_osp.use_log 

        if ret_levels is None:
            self.config_def['value'] = gas_osp.climatology_full_grid
            self.config_def['pressure'] = self.pressure_obj(gas_osp.pressure_full_grid)
        else:
            self.config_def['value'] = gas_osp.climatology_full_grid[ret_levels]
            self.config_def['pressure'] = self.pressure_obj(gas_osp.pressure_full_grid[ret_levels])

        return absorber.AbsorberVmrLevel.create(self, gas_name, **kwargs)

class CovarianceMUSES(CreatorMUSES, retrieval.CovarianceByComponent):

    num_channels = param.Scalar(int)

    def absorber_cov(self, gas_name, ret_type):
        gas_osp = OSP(gas_name, **self.osp_common())

        if ret_type == "log" and not gas_osp.use_log:
            raise Exception(f"OSP constraint for {gas_name} is for a log retrieval, but gas is not set up for a log retrieval")

        if gas_osp.constraint_filename is not None:
            return np.linalg.inv(gas_osp.constraint_matrix)
        else:
            raise param.ParamError(f"No constraint file found for gas: {gas_name}")

    def surface_temp_cov(self):
        tsur_osp = OSP("TSUR", **self.osp_common())
        return np.identity(self.num_channels()) * tsur_osp.covariance_matrix

    def surface_press_cov(self):
        psur_osp = OSP("PSUR", **self.osp_common())
        return psur_osp.covariance_matrix

    def temperature_cov(self):
        temp_osp = OSP("TATM", **self.osp_common())
        return np.linalg.inv(temp_osp.constraint_matrix)

    def create(self, **kwargs):

        # Get existing covariance values
        cov_values = self.values()

        # Add values based on retrieval type
        for rc_name in self.retrieval_components().keys():
            # Skip existing defined covariuance values
            if rc_name in cov_values:
                continue

            gas_match = re.match("absorber_levels/(.*)/(.*)", rc_name)
            if gas_match:
                ret_type = gas_match.groups()[0]
                gas_name = gas_match.groups()[1]

                cov_values[rc_name] = self.absorber_cov(gas_name, ret_type)

            elif rc_name == "surface_temperature":
                cov_values[rc_name] = self.surface_temp_cov()

            elif rc_name == "temperature_level":
                cov_values[rc_name] = self.temperature_cov()

            else:
                raise param.ParamError(f"Do not know how to set up covariance for {rc_name} retrieval component")

        self.config_def['values'] = dict(cov_values)

        return retrieval.CovarianceByComponent.create(self, **kwargs)
