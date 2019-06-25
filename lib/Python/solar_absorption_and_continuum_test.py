from .test_support import *
import os


# We are in the process of changing the solar model, so temporarily
# skip this test.
@skip
def test_solar_spectrum():
    spectrum_file = unit_test_data + "in/l1b/spec/spectra.dat"
    sounding_file = unit_test_data + "in/l1b/soundinginfo.dat"
    noise = rf.PrecomputedNoiseModel(HeritageFile(spectrum_file))
    l1b = rf.Level1bHeritage(sounding_file, spectrum_file, noise)
    h = rf.HdfFile(unit_test_data + "l2_fixed_level_static_input.h5")
    solar = rf.VectorSolarModel(h, l1b, "Solar")[0]
    wn = [12929.94, 12979.93, 13029.93, 13079.93, 13129.93,
          13179.93]
    expected = np.array([0.0726478736339, 0.0732139705623,
                         0.0727836689145, 0.073109528699,
                         0.0729386962257, 0.0727803937312])
    sd = rf.SpectralDomain(wn)
    spec = solar.solar_spectrum(sd)
    assert_array_almost_equal(spec.spectral_range().data(), expected)
    
