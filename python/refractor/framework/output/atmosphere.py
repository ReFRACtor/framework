import logging

import numpy as np

import refractor.framework as rf
from .base import OutputBase

logger = logging.getLogger(__name__)

class AtmosphereOutputBase(OutputBase):

    def __init__(self, output, step_index, atmosphere):
        super().__init__(output)

        self.step_index = step_index
        self.atm = atmosphere

    def output_atmosphere(self, iter_index):
        base_group_name = self.iter_step_group_name(self.step_index, iter_index)

        # Dimensions

        layer_dim = self.create_dimension("num_layer", self.atm.number_layer, self.step_index)
        level_dim = self.create_dimension("num_level", self.atm.number_layer+1, self.step_index)
        gas_dim = self.create_dimension("num_gas", self.atm.absorber.number_species, self.step_index)
        if self.atm.aerosol is not None:
            aer_dim = self.create_dimension("num_aerosol", self.atm.aerosol.number_particle, self.step_index)
        gas_str_dim = self.create_dimension("gas_string", step_index=self.step_index)
        aer_str_dim = self.create_dimension("aer_string", step_index=self.step_index)

        # Top level group

        atm_group = self.create_group(base_group_name + "Atmosphere")

        # Pressure

        pressure_var = self.create_variable("pressure_levels", atm_group, float, (level_dim,))
        pressure_var[:] = self.atm.pressure.pressure_grid().value.value
        pressure_var.units = self.atm.pressure.pressure_grid().units.name

        # Temperature
        temperature_var = self.create_variable("temperature_levels", atm_group, float, (level_dim,))
        temperature_var[:] = self.atm.temperature.temperature_grid(self.atm.pressure).value.value
        temperature_var.units = self.atm.temperature.temperature_grid(self.atm.pressure).units.name

        # Absorber
 
        absorber_group = self.create_group("Absorber", atm_group)

        gas_name_list = [ self.atm.absorber.gas_name(idx) for idx in range(self.atm.absorber.number_species) ]
        self.set_string_variable("gas_names", absorber_group, (gas_dim, gas_str_dim), gas_name_list)

        vmr_var = self.create_variable("vmr", absorber_group, float, (gas_dim, level_dim))

        # Which sensor to use for output of density information
        density_sensor_idx = 0

        gas_column = self.create_variable("gas_number_density_layer", absorber_group, float, (gas_dim, layer_dim))

        gas_density_layer = self.atm.absorber.gas_number_density_layer(density_sensor_idx)

        for gas_idx, gas_name in enumerate(gas_name_list):
            avmr = self.atm.absorber.absorber_vmr(gas_name)
            vmr_var[gas_idx, :] = avmr.vmr_grid(self.atm.pressure).value
            
            gas_column[gas_idx, :] = gas_density_layer.value.value[:, gas_idx]
            gas_column.units = gas_density_layer.units.name

        dry_air = self.create_variable("total_dry_air_number_density_layer", absorber_group, float, (layer_dim,))
        dry_air[:] = self.atm.absorber.total_air_number_density_layer(density_sensor_idx).value.value
        dry_air.units = self.atm.absorber.total_air_number_density_layer(density_sensor_idx).units.name
    
        # Aerosol 

        if self.atm.aerosol is not None:
            aerosol_group = self.create_group("Aerosol", atm_group)

            aer_name_list = self.atm.aerosol.aerosol_name
            self.set_string_variable("aerosol_names", aerosol_group, (aer_dim, aer_str_dim), aer_name_list)

            ext_var = self.create_variable("extinction", aerosol_group, float, (aer_dim, level_dim))

            for aer_idx, aer_name in enumerate(aer_name_list):
                aer_obj = self.atm.aerosol.aerosol_extinction(aer_idx)
                ext_var[aer_idx, :] = aer_obj.aerosol_extinction().value

class AtmosphereOutputRetrieval(rf.ObserverIterativeSolver, AtmosphereOutputBase):

    def __init__(self, output, step_index, atmosphere):
        rf.ObserverIterativeSolver.__init__(self)

        AtmosphereOutputBase.__init__(self, output, step_index, atmosphere)

    def notify_update(self, solver):
        try:
            iter_index = solver.num_accepted_steps
            self.output_atmosphere(iter_index)
        except Exception as exc:
            import traceback
            import sys
            exc_type, exc_value, exc_traceback = sys.exc_info()
            traceback.print_tb(exc_traceback)
            raise exc

class AtmosphereOutputSimulation(AtmosphereOutputBase):

    def __init__(self, output, step_index, atmosphere):
        super().__init__(output, step_index, atmosphere)

        self.output_atmosphere(None)
