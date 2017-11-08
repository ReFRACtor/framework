import os

import numpy as np

from .base import Creator, ParamPassThru
from .apriori import CreatorApriori
from .. import param

from refractor import framework as rf

class GasVmrAprioriMetL1b(Creator):
    "Creates a VMR apriori for a gas species using the TCCON method"


    l1b = param.InstanceOf(rf.Level1b)
    met = param.InstanceOf(rf.Meteorology)
    pressure = param.InstanceOf(rf.Pressure)
    altitudes = param.ObjectVector("altitude")
    reference_atm_file = param.Scalar(str)
    gas_name = param.Scalar(str)
    temp_avg_window = param.Scalar(int, default=11)

    def create(self, **kwargs):
        ref_atm_data = rf.HdfFile(self.reference_atm_file())
        apriori_obj = rf.GasVmrApriori(self.met(), self.l1b(), self.altitudes()[0], ref_atm_data, "/Reference_Atmosphere", self.gas_name(), self.temp_avg_window())
        return apriori_obj.apriori_vmr(self.pressure())

class GasVmrApriori(Creator):
    "Creates a VMR apriori for a gas species using the TCCON method"

    pressure = param.InstanceOf(rf.Pressure)
    temperature = param.InstanceOf(rf.Temperature)
    latitude = param.ArrayWithUnit(dims=1)
    time = param.Iterable(rf.Time)
    altitudes = param.ObjectVector("altitude")
    reference_atm_file = param.Scalar(str)
    gas_name = param.Scalar(str)
    temp_avg_window = param.Scalar(int, default=11)

    def create(self, **kwargs):
        ref_atm_data = rf.HdfFile(self.reference_atm_file())

        pressure_levels = self.pressure().pressure_grid.value.value
        temperature_levels = self.temperature().temperature_grid(self.pressure()).value.value

        apriori_obj = rf.GasVmrApriori(pressure_levels, temperature_levels, self.latitude().value[0], self.time()[0], \
                self.altitudes()[0], ref_atm_data, "/Reference_Atmosphere", self.gas_name(), self.temp_avg_window())
        return apriori_obj.apriori_vmr()

class AbsorberVmrLevel(CreatorApriori):
    "Creates a AbsorberVmrLevel that supplies a AbsorberVmr class for use in an creating an Atmosphere"

    pressure = param.InstanceOf(rf.Pressure)

    def create(self, gas_name=None, **kwargs):

        if gas_name is None:
            raise param.ParamError("gas_name not supplied to creator %s" % self.__class__.__name__)

        return rf.AbsorberVmrLevel(self.pressure(), self.apriori(), self.retrieval_flag(), gas_name)

class AbsorberVmrMet(CreatorApriori):

    met = param.InstanceOf(rf.Meteorology)
    pressure = param.InstanceOf(rf.Pressure)

    def create(self, gas_name=None, **kwargs):

        if gas_name is None:
            raise param.ParamError("gas_name not supplied to creator %s" % self.__class__.__name__)

        return rf.AbsorberVmrMet(self.met(), self.pressure(), self.apriori()[0], bool(self.retrieval_flag()[0]), gas_name)

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
    altitudes = param.ObjectVector("altitude")
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

