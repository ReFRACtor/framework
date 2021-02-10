import logging

import numpy as np
import netCDF4

from refractor import framework as rf
from .base import OutputBase

logger = logging.getLogger(__name__)

class ForwardModelRadianceOutput(rf.ObserverPtrNamedSpectrum, OutputBase):

    # Maps the NamedSpectrum names to groups where to save data
    spectrum_type_group_names = {
        "high_res_rt": "Monochromatic",
        "convolved": "Convolved",
    }

    dimension_names = {
        "high_res_rt": "monochromatic_grid_s{}c{}",
        "convolved": "instrument_grid_s{}c{}",
    }

    def __init__(self, output, step_index, solver=None):
        # Required to initialize director
        rf.ObserverPtrNamedSpectrum.__init__(self)

        self.output = output
        self.step_index = step_index

        self.solver = solver

    def notify_update(self, named_spectrum):
        if named_spectrum.name not in self.spectrum_type_group_names:
            logger.debug("Ignoring non saved spectrum type: {}".format(named_spectrum.name))
            return

        # Determine where to store values
        if self.solver is not None:
            # Add 1 because solver will not have incremented its internal number yet
            iter_index = self.solver.num_accepted_steps + 1
        else:
            iter_index = None

        fill_value = netCDF4.default_fillvals['f8']

        group_name = self.iter_step_group_name(self.step_index, iter_index)
        group_name += "Spectrum"
        group_name += "/" + "Channel_{}".format(named_spectrum.index + 1)
        group_name += "/" + self.spectrum_type_group_names[named_spectrum.name]

        logger.debug("Storing spectrum of type {} into output file at {}".format(named_spectrum.name, group_name))

        group = self.output.createGroup(group_name)

        dim_name = self.dimension_names[named_spectrum.name].format(self.step_index+1, named_spectrum.index+1)
        if not dim_name in self.output.dimensions:
            # Make radiance dimension unitless because the size could change between iterations due
            # to sample_grid retrieval
            self.output.createDimension(dim_name, None)

        grid_name = "grid"
        if grid_name in group.variables:
            grid_data = group[grid_name]

            # Populate with fill value for the case where data size changes when overwriting
            grid_data[:] = fill_value
        else:
            grid_data = group.createVariable(grid_name, float, (dim_name,), fill_value=fill_value)

        npoints = named_spectrum.spectral_domain.data.shape[0]
        grid_data[:npoints] = named_spectrum.spectral_domain.data

        radiance_name = "radiance"
        if radiance_name in group.variables:
            radiance_data = group[radiance_name]

            # Populate with fill value for the case where data size changes when overwriting
            radiance_data[:] = fill_value
        else:
            radiance_data = group.createVariable(radiance_name, float, (dim_name,), fill_value=fill_value)

        npoints = named_spectrum.spectral_range.data.shape[0]
        radiance_data[:npoints] = named_spectrum.spectral_range.data
        radiance_data.units = named_spectrum.spectral_range.units.name

        rad_ad = named_spectrum.spectral_range.data_ad
        if not rad_ad.is_constant:

            jac_dim = "jacobian_s{}".format(self.step_index+1)
            if jac_dim not in self.output.dimensions:
                self.output.createDimension(jac_dim, rad_ad.jacobian.shape[1])

            jacobian_name = "jacobian"
            if jacobian_name in group.variables:
                jacobian_data = group[jacobian_name]

                jacobian_data[:] = fill_value
            else:
                jacobian_data = group.createVariable(jacobian_name, float, (dim_name, jac_dim), fill_value=fill_value)

            jacobian_data[:] = rad_ad.jacobian[:]

class ObservationRadianceOutput(rf.ObserverIterativeSolver, OutputBase):
    "Outputs Level1B related values into the output file"

    def __init__(self, output, step_index, l1b, forward_model):
        # Required to initialize director
        rf.ObserverIterativeSolver.__init__(self)

        self.output = output
        self.step_index = step_index

        self.l1b = l1b
        self.forward_model = forward_model

        self.write_observation()

    def write_observation(self):

        # Don't write the observation values twice
        if "Observation" in self.output.groups:
            return

        obs_group = self.output.createGroup("Observation")

        for channel_idx in range(self.l1b.number_spectrometer()):
            channel_group = obs_group.createGroup("Channel_{}".format(channel_idx+1))

            channel_grid = self.l1b.sample_grid(channel_idx).data
            channel_rad  = self.l1b.radiance(channel_idx).data
            channel_unc  = self.l1b.radiance(channel_idx).uncertainty

            obs_dim = "observation_grid_s{}c{}".format(self.step_index+1, channel_idx+1)
            if not obs_dim in self.output.dimensions:
                self.output.createDimension(obs_dim, channel_grid.shape[0])

            grid_units = self.l1b.sample_grid(channel_idx).units.name
            rad_units = self.l1b.radiance(channel_idx).units.name

            grid_name = "grid"
            obs_grid = channel_group.createVariable(grid_name, float, (obs_dim,))
            obs_grid[:] = channel_grid
            obs_grid.units = grid_units

            rad_name = "radiance"
            obs_rad = channel_group.createVariable(rad_name, float, (obs_dim,))
            obs_rad[:] = channel_rad
            obs_rad.units = rad_units

            unc_name = "uncertainty"
            obs_unc = channel_group.createVariable(unc_name, float, (obs_dim,))
            obs_unc[:] = channel_unc

    def write_sample_list(self, spectrum_group):
        "Update sample list for each iteration/step"

        for channel_idx in range(self.l1b.number_spectrometer()):
            channel_group = spectrum_group.createGroup("Channel_{}".format(channel_idx+1))

            channel_samples = np.array(self.forward_model.spectral_grid.pixel_list(channel_idx))

            if channel_samples.shape[0] > 0:
                inst_dim = "instrument_grid_s{}c{}".format(self.step_index+1, channel_idx+1)
                if not inst_dim in self.output.dimensions:
                    self.output.createDimension(inst_dim, channel_samples.shape[0])

                if "observation_sample_indexes" in channel_group.variables:
                    obs_indexes = channel_group["observation_sample_indexes"]
                else:
                    obs_indexes = channel_group.createVariable("observation_sample_indexes", int, (inst_dim,))

                obs_indexes[:] = channel_samples

    def notify_update(self, solver):
        iter_index = solver.num_accepted_steps
        base_group_name = self.iter_step_group_name(self.step_index, iter_index)
        
        group_name = base_group_name + "Spectrum"

        group = self.output.createGroup(group_name)

        self.write_sample_list(group)
