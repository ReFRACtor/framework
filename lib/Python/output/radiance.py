import logging

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

        self.create_dimensions()

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

        if "grid" in group.variables:
            grid = group['grid']
        else:
            grid = group.createVariable("grid", float, (self.dimension_names[named_spectrum.name],))
        grid[:] = named_spectrum.spectral_domain.data

        if "radiance" in group.variables:
            radiance = group['radiance']
        else:
            radiance = group.createVariable("radiance", float, (self.dimension_names[named_spectrum.name],))
        radiance[:] = named_spectrum.spectral_range.data
