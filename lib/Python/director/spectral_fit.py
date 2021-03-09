import numpy as np

from refractor import framework as rf

class SpectralFitSampleGrid(rf.SampleGridImpBase):
    """Applies a retrievable shift and squeeze to the instrument's sample grid."""
    
    def __init__(self, shift : float, squeeze : float, instrument_grid : rf.SpectralDomain, sensor_name : str):
        retrieval_coeffs = [shift, squeeze]
        rf.SampleGridImpBase.__init__(self, retrieval_coeffs)
        self.instrument_grid = instrument_grid
        self.sensor_name = sensor_name

    # Since sample_grid is %python_attribute in rf.SampleGrid, we must implement the 
    # renamed function instead of the Python property otherwise an infinite loop will occur
    def _v_sample_grid(self) -> rf.SpectralDomain: 

        # Obtain updated coefficeints from retrieval vector
        mapped_coeffs = self.state_mapping.mapped_state(self.coefficient)

        # Combine mapped_coeffs with self.instrtument_grid in some way
        # implement shift and squeeze

        retrieved_grid = rf.ArrayAd_double_1(self.instrument_grid.data.shape[0], mapped_coeffs.number_variable)

        # Compute at each wavelength to ensure correct auto derivative handling
        for wvl_idx in range(retrieved_grid.rows):
            retrieved_grid[wvl_idx] = mapped_coeffs[0] * (1 + mapped_coeffs[1]) * self.instrument_grid.data[wvl_idx]

        # We need to add sample index or else the code things the measure radiance is empty
        sample_index = np.arange(retrieved_grid.rows)

        return rf.SpectralDomain(retrieved_grid, sample_index, self.instrument_grid.units)

    def _v_pixel_grid(self) -> rf.SpectralDomain:
        return self._v_sample_grid()

    def _v_sub_state_identifier(self) -> str:
        "Defines the name of the group of items this class represents in the retrieval vector, one per sensor"

        return f"spectral_fit/{self.sensor_name}"

    def state_vector_name_i(self, i):
        "Defines the names of the coefficients stored in the retrieval vector"

        coeff_names = [ 
           f"Instrument Grid Spectral Fit {self.sensor_name} Shift",
           f"Instrument Grid Spectral Fit {self.sensor_name} Squeeze"
        ]

        return coeff_names[i]

    def clone(self) -> rf.SampleGridImpBase:
        "Creates a new copy of this class"

        return self.__class__(self.coefficient[0].value, self.coefficient[1].value, self.instrument_grid)
