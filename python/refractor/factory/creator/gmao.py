from functools import lru_cache

from .base import Creator
from .. import param

from refractor import framework as rf

from refractor.input.gmao import GMAOSupplier

# Cache supplier creation for a given lon/lat/time to avoid unnecessary searching
@lru_cache()
def _create_supplier(longitude, latitude, when, base_dir):
    return GMAOSupplier(longitude, latitude, when, base_dir)

class CreatorGMAOBase(Creator):

    l1b = param.InstanceOf(rf.Level1b) 
    gmao_directory = param.Scalar(str)

    def gmao_supplier(self):

        l1b_obj = self.l1b()

        longitude = l1b_obj.longitude(0).value
        latitude = l1b_obj.latitude(0).value
        obs_time = l1b_obj.time(0).as_datetime()

        return _create_supplier(longitude, latitude, obs_time, self.gmao_directory())

class PressureGMAO(CreatorGMAOBase):

    def create(self, **kwargs):

        supplier = self.gmao_supplier()

        return rf.PressureSigma(supplier.pressure, supplier.surface_pressure)

class TemperatureProfileGMAO(CreatorGMAOBase):

    def create(self, **kwargs):

        supplier = self.gmao_supplier()

        return supplier.temperature
