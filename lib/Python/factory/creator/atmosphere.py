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
    #relative_humidty = 
    #ground = 
    #constant = 

    def create(self, **kwargs):

        pressure = self.common_store["pressure"] = self.pressure()
        temperature = self.common_store["temperature"] = self.temperature()
        altitudes = self.common_store["altitudes"] = self.altitudes()
        absorber = self.common_store["absorber"] = self.absorber()

        print("-->", absorber)

        #return fp.AtmosphereOco(absorber, pressure, temperature, 
        #                      rel_hum, ground, altitudes, fp.DefaultConstant())

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
    
    latitude = param.DoubleWithUnit()
    surface_height = param.DoubleWithUnit()

    pressure = param.InstanceOf(rf.Pressure)
    temperature = param.InstanceOf(rf.Temperature)

    num_channels = param.Scalar(int)

    def create(self, **kwargs):
        altitudes = rf.vector_altitude()
        for chan_idx in range(self.num_channels()):
            chan_alt = rf.AltitudeHydrostatic(self.pressure(), self.temperature(), self.latitude(channel_index=chan_idx), self.surface_height(channel_index=chan_idx))
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
    table_scale = param.Array(dims=1)

    def create(self):

        
        pass

class AbsorberGasDefinition(ParamPassThru):
    "Defines the interface expected for VMR config defnition blocks, values are pass through as a dictionary"

    vmr = param.InstanceOf(rf.AbsorberVmr)
    absco = param.InstanceOf(rf.Absco)

class AbsorberAbsco(Creator):
    "Creates an AbsorberAbsco object that statisfies the AtmosphereCreato;rs absorber value"

    gases = param.Iterable()
    pressure = param.InstanceOf(rf.Pressure)
    temperature = param.InstanceOf(rf.Temperature)
    altitudes = param.AnyValue() #param.ObjectVector() ### TODO
    num_sub_layers = param.Scalar(int, required=False)
 
    def create(self, **kwargs):

        vmrs = rf.vector_absorber_vmr
        abscos = rf.vector_absco

        for gas_name in self.gases():
            self.register_parameter(gas_name, param.Dict())
            gas_def = self.param(gas_name, gas_name=gas_name)
            
            vmrs.push_back(gas_def['vmr'])
            absco.push_back(gas_def['absco'])

        if self.num_sub_layers() is not None:
            return AbsorberAbsco(vmrs, self.pressure(), self.temperature(), self.altitudes(), abscos, self.constants(), self.number_sub_layers())
        else:
            return AbsorberAbsco(vmrs, self.pressure(), self.temperature(), self.altitudes(), abscos, self.constants())


class GasVMRFromConstant(Creator):
    def create(self, **kwargs):
        pass
