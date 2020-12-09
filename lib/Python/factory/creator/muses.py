import re
import logging

import numpy as np

from .base import Creator
from .. import param

from refractor import framework as rf
from . import atmosphere
from . import absorber
from . import retrieval

from refractor.muses_osp import PressureOSP, SpeciesOSP

logger = logging.getLogger(__name__)

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

    def pressure_obj(self, pressure_levels, surface_pressure=None):
        if surface_pressure is None:
            return rf.PressureSigma(pressure_levels, pressure_levels[-1])
        else:
            # Keep PressureSigma for remapping the levels 
            pressure_levels[-1] = surface_pressure
            return rf.PressureSigma(pressure_levels, surface_pressure)

class PressureGridMUSES(CreatorMUSES):

    pressure_levels = param.Choice(param.Array(dims=1), param.NoneValue())
    surface_pressure = param.Choice(param.Scalar(float), param.NoneValue())

    def create(self, **kwargs):

        pressure_levels = self.pressure_levels()
        surface_pressure = self.surface_pressure()

        if pressure_levels is None:
            # We do want a cut-off of 1000mb here for the purposes of creating
            # the sigma levels, this happens with the default keyword value for pressure_cutoff
            press_osp = PressureOSP(**self.osp_common())

            logger.debug(f"Loading pressure grid from: {press_osp.fm_pressure_filename}")

            pressure_levels = press_osp.fm_pressure_grid

        if surface_pressure is None:
            # Disable cutoff to get all pressure levels present
            press_osp = PressureOSP(pressure_cutoff=None, **self.osp_common())

            fm_press = press_osp.fm_pressure_grid
            press_obj = rf.PressureSigma(fm_press, fm_press[-1])

            temp_osp = SpeciesOSP("TATM", pressure_cutoff=None, **self.osp_common())
            fm_temp = temp_osp.fm_climatology
            
            surf_press_def = {
                'pressure': press_obj,
                'temperature_profile': fm_temp,
                'latitude': self.l1b().latitude(0),
                'surface_height': self.l1b().altitude(0),
            }

            logger.debug("Determining surface pressure from altitude profile")

            surf_press_creator = atmosphere.SurfacePressureFromAltitude(surf_press_def)
            surface_pressure = surf_press_creator.create()

        # Ensure sigma levels does not remap pressure levels by explicitly setting the
        # surface pressure here
        pressure_levels[-1] = surface_pressure

        press_obj = rf.PressureSigma(pressure_levels, surface_pressure)

        logger.debug(f"Using forward model pressure grid:\n{press_obj.pressure_grid.value.value}")

        return press_obj

class TemperatureMUSES(CreatorMUSES, atmosphere.TemperatureLevel):

    temperature_profile = param.Choice(atmosphere.TemperatureLevel.temperature_profile, param.NoneValue())

    def create(self, **kwargs):

        tatm_osp = SpeciesOSP("TATM", **self.osp_common())
        ret_levels = tatm_osp.retrieval_levels
        osp_pressure_grid = tatm_osp.fm_pressure_grid[ret_levels]

        if self.temperature_profile() is None:

            logger.debug(f"Loading temperature profile from: {tatm_osp.climatology_filename}")

            self.config_def["temperature_profile"] = tatm_osp.fm_climatology[ret_levels]

        else:

            logger.debug(f"Resampling supplied pressure profile to pressure grid: {osp_pressure_grid}")

            # Resample supplied value so that covariance from MUSES can be used
            inp_pressure_grid = self.pressure().pressure_grid.value.value

            self.config_def["temperature_profile"] = np.interp(osp_pressure_grid, inp_pressure_grid, self.temperature_profile())

        # Either way the pressure grid used will be the OSP pressure grid
        self.config_def["pressure"] = self.pressure_obj(osp_pressure_grid)

        return atmosphere.TemperatureLevel.create(self, **kwargs)

class SurfaceTemperatureMUSES(CreatorMUSES, atmosphere.SurfaceTemperature):

    num_channels = param.Scalar(int)
    surface_temperature = param.Choice(atmosphere.SurfaceTemperature.surface_temperature, param.NoneValue())

    def create(self, **kwargs):

        tsur_osp = SpeciesOSP("TSUR", **self.osp_common())

        logger.debug(f"Loading surface temperature from: {tsur_osp.climatology_filename}")

        if self.surface_temperature() is None:
            self.config_def['surface_temperature'] = rf.ArrayWithUnit_double_1(np.full(self.num_channels(), tsur_osp.species_climatology), "K")

        return atmosphere.SurfaceTemperature.create(self, **kwargs)

class AbsorberVmrMUSES(CreatorMUSES, absorber.AbsorberVmrLevel):

    vmr_profile = param.Choice(absorber.AbsorberVmrLevel.vmr_profile, param.NoneValue(), required=False)
    mapping = param.Choice(absorber.AbsorberVmrLevel.mapping, param.NoneValue(), required=False)

    def create(self, gas_name, **kwargs):
        gas_osp = SpeciesOSP(gas_name, **self.osp_common())

        if self.vmr_profile() is None:
            ret_levels = gas_osp.retrieval_levels

            if ret_levels is not None:
                logger.debug(f"Loading {len(ret_levels)} retrieval levels for {gas_name} from {gas_osp.climatology_filename}")
                self.config_def['vmr_profile'] = gas_osp.fm_climatology[ret_levels]
            else:
                logger.debug(f"Loading full forward model levels for {gas_name} from {gas_osp.climatology_filename}")
                self.config_def['vmr_profile'] = gas_osp.fm_climatology

        if self.mapping() is None:
            surface_pressure = self.pressure().surface_pressure.value.value

            if ret_levels is not None:
                pressure_from = self.pressure_obj(gas_osp.fm_pressure_grid[ret_levels])
            else:
                pressure_from = self.pressure_obj(gas_osp.fm_pressure_grid)

            if gas_osp.use_log:
                mapping_first = rf.StateMappingLog()
            else:
                mapping_first = rf.StateMappingLinear()

            logger.debug(f"Mapping {gas_name} levels from retrieval pressure grid:\n{pressure_from.pressure_grid.value.value}\nto forward model pressure grid:\n{self.pressure().pressure_grid.value.value}")

            mapping_interp = rf.StateMappingInterpolateLogLog(self.pressure(), pressure_from)

            mappings = rf.vector_state_mapping()
            mappings.push_back(mapping_first)
            mappings.push_back(mapping_interp)

            mapping_composite = rf.StateMappingComposite(mappings)

            self.config_def['mapping'] = mapping_composite
            self.config_def['coeff_pressure'] = pressure_from

        return absorber.AbsorberVmrLevel.create(self, gas_name, **kwargs)

class CovarianceMUSES(CreatorMUSES, retrieval.CovarianceByComponent):

    num_channels = param.Scalar(int)

    def covariance(self, constraint):

        return np.linalg.inv(constraint)

    def absorber_cov(self, gas_name, ret_type):
        gas_osp = SpeciesOSP(gas_name, **self.osp_common())

        if ret_type == "log" and not gas_osp.use_log:
            raise Exception(f"OSP constraint for {gas_name} is for a log retrieval, but gas is not set up for a log retrieval")

        if gas_osp.constraint_filename is not None:
            logger.debug(f"Loading constraint for {gas_name} from {gas_osp.constraint_filename}")
            return self.covariance(gas_osp.constraint_matrix)
        else:
            raise param.ParamError(f"No constraint file found for gas: {gas_name}")

    def surface_temp_cov(self):
        tsur_osp = SpeciesOSP("TSUR", **self.osp_common())
        logger.debug(f"Loading constraint for surface temperature from {tsur_osp.constraint_filename}")
        return np.identity(self.num_channels()) * tsur_osp.covariance_matrix

    def surface_press_cov(self):
        psur_osp = SpeciesOSP("PSUR", **self.osp_common())
        logger.debug(f"Loading constraint for surface pressure from {psur_osp.constraint_filename}")
        return psur_osp.covariance_matrix

    def temperature_cov(self):
        temp_osp = SpeciesOSP("TATM", **self.osp_common())
        logger.debug(f"Loading constraint for temperature from {temp_osp.constraint_filename}")
        return self.covariance(temp_osp.constraint_matrix)

    def create(self, **kwargs):

        # Get existing covariance values
        storage = self.interstep_storage()
        if storage is not None:
            if len(storage) == 0:
                storage.update(self.values())
            cov_values = storage
        else:
            cov_values = self.values()

        # Add values based on retrieval type
        for rc_name in self.retrieval_components().keys():
            # Skip existing defined covariance values
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
