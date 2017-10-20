import os

import numpy as np

from .base import Creator, ParamPassThru
from .apriori import CreatorApriori
from .. import param

from refractor import framework as rf

class AtmosphereCreator(Creator):
    "Creates an atmosphere object"

    pressure = param.InstanceOf(rf.Pressure)
    temperature = param.InstanceOf(rf.Temperature)
    altitudes = param.AnyValue() #param.ObjectVector() ### TODO
    absorber = param.InstanceOf(rf.Absorber)
    relative_humidity = param.InstanceOf(rf.RelativeHumidity)
    ground = param.InstanceOf(rf.Ground)
    #aerosol =
    constants = param.InstanceOf(rf.Constant)

    def create(self, **kwargs):

        pressure = self.common_store["pressure"] = self.pressure()
        temperature = self.common_store["temperature"] = self.temperature()
        altitudes = self.common_store["altitudes"] = self.altitudes()
        absorber = self.common_store["absorber"] = self.absorber()
        relative_humidity = self.relative_humidity()
        ground = self.ground()

        return rf.AtmosphereOco(absorber, pressure, temperature, relative_humidity, ground, altitudes, self.constants())

class PressureSigma(CreatorApriori):
    "Creates a PressureSigma object statisfying the AtmosphereCreator's pressure parameter"

    a_coeff = param.Array(dims=1)
    b_coeff = param.Array(dims=1)

    def create(self, **kwargs):
        # ap and flag loaded as arrays, so just get first value
        ap_psurf = self.apriori()[0]
        ret_flag = bool(self.retrieval_flag()[0])

        return rf.PressureSigma(self.a_coeff(), self.b_coeff(), ap_psurf, ret_flag)

class TemperatureMet(CreatorApriori):
    "Creates a TemperatureMet object statisfying the AtmosphereCreator's temperature parameter"
    
    met = param.InstanceOf(rf.Meteorology)
    pressure = param.InstanceOf(rf.Pressure)

    def create(self, **kwargs):
        ap_offset = self.apriori()[0]
        ret_flag = bool(self.retrieval_flag()[0])

        return rf.TemperatureMet(self.met(), self.pressure(), ap_offset, ret_flag)

class AltitudeHydrostatic(Creator):
    "Creates a AltitudeHydrostatic object statisfying the AtmosphereCreator's altitude parameter"
    
    latitude = param.ArrayWithUnit(dims=1)
    surface_height = param.ArrayWithUnit(dims=1)

    pressure = param.InstanceOf(rf.Pressure)
    temperature = param.InstanceOf(rf.Temperature)

    num_channels = param.Scalar(int)

    def create(self, **kwargs):
        # These are per channel
        latitudes = self.latitude()
        surface_heights = self.surface_height()

        altitudes = rf.vector_altitude()
        for chan_idx in range(self.num_channels()):
            chan_alt = rf.AltitudeHydrostatic(self.pressure(), self.temperature(), latitudes[chan_idx], surface_heights[chan_idx])
            altitudes.push_back(chan_alt)

        return altitudes

class GasVmrApriori(Creator):
    "Creates a VMR apriori for a gas species using the TCCON method"


    l1b = param.InstanceOf(rf.Level1b)
    met = param.InstanceOf(rf.Meteorology)
    pressure = param.InstanceOf(rf.Pressure)
    altitudes = param.AnyValue() #param.ObjectVector() ### TODO
    reference_atm_file = param.Scalar(str)
    gas_name = param.Scalar(str)
    temp_avg_window = param.Scalar(int, default=11)

    def create(self, **kwargs):
        ref_atm_data = rf.HdfFile(self.reference_atm_file())
        apriori_obj = rf.GasVmrApriori(self.met(), self.l1b(), self.altitudes()[0], ref_atm_data, "/Reference_Atmosphere", self.gas_name(), self.temp_avg_window())
        return apriori_obj.apriori_vmr(self.pressure())

class AbsorberVmrLevel(CreatorApriori):
    "Creates a AbsorberVmrLevel that supplies a AbsorberVmr class for use in an creating an Atmosphere"

    pressure = param.InstanceOf(rf.Pressure)

    def create(self, gas_name=None, **kwargs):

        if gas_name is None:
            raise param.ParamError("gas_name not supplied to creator %s" % self.__class__.__name__)

        return rf.AbsorberVmrLevel(self.pressure(), self.apriori(), self.retrieval_flag(), gas_name)

class AbscoHdf(Creator):

    absco_base_path = param.Scalar(str)
    filename = param.Scalar(str)
    table_scale = param.Choice(param.Iterable(), param.Scalar(float), default=1.0)
    spec_win = param.InstanceOf(rf.SpectralWindow)

    def create(self, **kwargs):

        absco_filename = os.path.join(self.absco_base_path(), self.filename())

        if not os.path.exists(absco_filename):
            raise param.ParamError("HDF ABSCO filename does not exist: %s" % absco_filename)

        table_scale = self.table_scale()

        if np.isscalar(table_scale):
            return rf.AbscoHdf(absco_filename, table_scale)
        else:
            spectral_bound = self.spec_win().spectral_bound

            # Convert to vector to match interface
            table_scale_vector = rf.vector_double()
            for val in table_scale:
                table_scale_vector.push_back(val)

            return rf.AbscoHdf(absco_filename, spectral_bound, table_scale_vector)

class AbsorberGasDefinition(ParamPassThru):
    "Defines the interface expected for VMR config defnition blocks, values are pass through as a dictionary"

    vmr = param.InstanceOf(rf.AbsorberVmr)
    absorption = param.InstanceOf(rf.GasAbsorption)

class AbsorberAbsco(Creator):
    "Creates an AbsorberAbsco object that statisfies the AtmosphereCreato;rs absorber value"

    gases = param.Iterable()
    pressure = param.InstanceOf(rf.Pressure)
    temperature = param.InstanceOf(rf.Temperature)
    altitudes = param.AnyValue() #param.ObjectVector() ### TODO
    num_sub_layers = param.Scalar(int, required=False)
    constants = param.InstanceOf(rf.Constant)
 
    def create(self, **kwargs):

        vmrs = rf.vector_absorber_vmr()
        absorptions = rf.vector_gas_absorption()

        for gas_name in self.gases():
            self.register_parameter(gas_name, param.Dict())
            gas_def = self.param(gas_name, gas_name=gas_name)

            if not "vmr" in gas_def:
                raise param.ParamError("vmr value not in gas definition for gas: %s" % gas_name)

            if not "absorption" in gas_def:
                raise param.ParamError("absorption value not in gas definition for gas: %s" % gas_name)

            vmrs.push_back(gas_def['vmr'])
            absorptions.push_back(gas_def['absorption'])

        if self.num_sub_layers() is not None:
            return rf.AbsorberAbsco(vmrs, self.pressure(), self.temperature(), self.altitudes(), absorptions, self.constants(), self.number_sub_layers())
        else:
            return rf.AbsorberAbsco(vmrs, self.pressure(), self.temperature(), self.altitudes(), absorptions, self.constants())


class GasVMRFromConstant(Creator):
    def create(self, **kwargs):
        pass

class RelativeHumidity(Creator):

    pressure = param.InstanceOf(rf.Pressure)
    temperature = param.InstanceOf(rf.Temperature)
    absorber = param.InstanceOf(rf.Absorber)

    def create(self, **kwargs):
        return rf.RelativeHumidity(self.absorber(), self.temperature(), self.pressure())
