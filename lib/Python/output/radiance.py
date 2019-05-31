import logging

from refractor import framework as rf

logger = logging.getLogger(__name__)

class ForwardModelRadianceOutput(rf.ObserverPtrNamedSpectrum):

    # Maps the NamedSpectrum names to groups where to save data
    spectrum_type_group_names = {
        "high_res_rt": "Monochromatic",
        "convolved": "Convolved",
    }

    spectrum_type_dimension_names = {
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

    def create_dimensions(self):

        self.dimensions = {}
        for spec_type, dim_name in self.spectrum_type_dimension_names.items():
            if dim_name in self.output.dimensions:
                self.dimensions[spec_type] = self.output.dimensions[dim_name]
            else:
                self.dimensions[spec_type] = self.output.createDimension(dim_name, 0)

    def notify_update(self, named_spectrum):

        if named_spectrum.name not in self.spectrum_type_group_names:
            logger.debug("Ignoring non saved spectrum type: {}".format(named_spectrum.name))
            return

        # Determine where to store values
        group_name = "Spectrum"
        if self.step_index is not None:
            group_name += "/Step_{}".format(self.step_index)

        # See if we can determine and iteration index for retrievals
        if self.solver is not None:
            # Add 1 because solver will not have incremented its internal number yet
            iter_index = self.solver.num_accepted_steps + 1
            if iter_index == 0:
                group_name += "/Initial"
            else:
                group_name += "/Iteration_{}".format(iter_index)

        group_name += "/" + self.spectrum_type_group_names[named_spectrum.name]

        logger.debug("Storing spectrum of type {} into output file at {}".format(named_spectrum.name, group_name))

        group = self.output.createGroup(group_name)

        if "grid" in group.variables:
            grid = group['grid']
        else:
            grid = group.createVariable("grid", float, (self.spectrum_type_dimension_names[named_spectrum.name],))
        grid[:] = named_spectrum.spectral_domain.data

        if "radiance" in group.variables:
            radiance = group['radiance']
        else:
            radiance = group.createVariable("radiance", float, (self.spectrum_type_dimension_names[named_spectrum.name],))
        radiance[:] = named_spectrum.spectral_range.data
