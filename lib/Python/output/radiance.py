import logging

import numpy as np

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
        "high_res_rt": "monochromatic_grid",
        "convolved": "instrument_grid",
    }

    def __init__(self, output, step_index, solver=None):
        # Required to initialize director
        rf.ObserverPtrNamedSpectrum.__init__(self)

        self.output = output
        self.step_index = step_index

        self.solver = solver

        self.create_dimensions(self.step_index)

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

        group_name = self.iter_step_group_name(self.step_index, iter_index)
        group_name += "Spectrum"
        group_name += "/" + self.spectrum_type_group_names[named_spectrum.name]

        logger.debug("Storing spectrum of type {} into output file at {}".format(named_spectrum.name, group_name))

        group = self.output.createGroup(group_name)

        dim_name = self.dimensions[named_spectrum.name].name

        if "grid" in group.variables:
            grid = group['grid']
        else:
            grid = group.createVariable("grid", float, (dim_name,))
        grid[:] = named_spectrum.spectral_domain.data
        grid.units = named_spectrum.spectral_domain.units.name

        if "radiance" in group.variables:
            radiance = group['radiance']
        else:
            radiance = group.createVariable("radiance", float, (dim_name,))
        radiance[:] = named_spectrum.spectral_range.data
        radiance.units = named_spectrum.spectral_range.units.name

class ObservationRadianceOutput(rf.ObserverIterativeSolver, OutputBase):
    "Outputs Level1B related values into the output file"

    dimension_names = {
        "observation_grid": "observation_grid",
        "instrument_grid": "instrument_grid",
    }

    def __init__(self, output, step_index, l1b, forward_model):
        # Required to initialize director
        rf.ObserverIterativeSolver.__init__(self)

        self.output = output
        self.step_index = step_index

        self.l1b = l1b
        self.forward_model = forward_model

        self.create_dimensions(self.step_index)

        self.write_observation()

    def write_observation(self):

        # Don't write the observation values twice
        if "Observation" in self.output.groups:
            return

        channel_grid = []
        channel_rad = []
        channel_unc = []

        grid_units = None
        rad_units = None
        for channel_idx in range(self.l1b.number_spectrometer()):
            channel_grid.append(self.l1b.sample_grid(channel_idx).data)
            channel_rad.append(self.l1b.radiance(channel_idx).data)
            channel_unc.append(self.l1b.radiance(channel_idx).uncertainty)

            grid_units = self.l1b.sample_grid(channel_idx).units.name
            rad_units = self.l1b.radiance(channel_idx).units.name

        group = self.output.createGroup("Observation")

        obs_dim = self.dimensions["observation_grid"].name

        obs_grid = group.createVariable("grid", float, (obs_dim,))
        obs_grid[:] = np.concatenate(channel_grid)
        obs_grid.units = grid_units

        obs_rad = group.createVariable("radiance", float, (obs_dim,))
        obs_rad[:] = np.concatenate(channel_rad)
        obs_rad.units = rad_units

        obs_unc = group.createVariable("uncertainty", float, (obs_dim,))
        obs_unc[:] = np.concatenate(channel_unc)

    def write_sample_list(self, group):
        "Update sample list for each iteration/step"

        channel_samples = []
        for channel_idx in range(self.l1b.number_spectrometer()):
            channel_samples.append(np.array(self.forward_model.spectral_grid.pixel_list(channel_idx)))

        if "observation_sample_indexes" in group.variables:
            obs_indexes = group["observation_sample_indexes"]
        else:
            inst_dim = self.dimensions["instrument_grid"].name
            obs_indexes = group.createVariable("observation_sample_indexes", int, (inst_dim,))

        obs_indexes[:] = np.concatenate(channel_samples)

    def notify_update(self, solver):

        iter_index = solver.num_accepted_steps
        base_group_name = self.iter_step_group_name(self.step_index, iter_index)
        
        group_name = base_group_name + "Spectrum"

        group = self.output.createGroup(group_name)

        self.write_sample_list(group)

