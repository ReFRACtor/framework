from test_support import *
from refractor_swig import *
import numpy as np

def test_cxx_serialization():
    '''Base test with ILS serialization, where class is all C++'''
    # These values come from the IlsGrating C++ unit test, just to give us
    # something reasonable to work with.
    d = DispersionPolynomial([1.28695614e+04, 1.99492886e-01],
                             Unit("cm^-1"),
                             np.arange(1,1806,1).astype(np.double),
                             "Test band")
    hf = HdfFile(unit_test_data + "in/ils/ils_linear_table.h5")
    ils_func = IlsTableLinear(hf, 0, "A-Band", "o2")
    ils = IlsGrating(d, ils_func)
    t = serialize_write_string(ils)
    ilsr = serialize_read_generic_string(t)
    # Not really a test that these are identical, but the actual serialization
    # gets tested at the C++ level. We just want to make sure that the
    # IlsImpBase stuff works with C++ derived classes and python director
    # derived ones.
    assert str(ils) == str(ilsr)

class IlsInPython(IlsImpBase):
    '''A dummy class to test that we serialize this correctly'''
    def __init__(self):
        self.data = 10
        d = DispersionPolynomial([1.28695614e+04, 1.99492886e-01],
                                 Unit("cm^-1"),
                                 np.arange(1,1806,1).astype(np.double),
                                 "Test band")
        super().__init__(d, DoubleWithUnit(0.0, "cm^-1"))

    def apply_ils(self, wn, rad, pixel_list):
        print("Dummy ils, with data %d" % self.data)
        return np.array([1.0,2.0,3.0])

# Doesn't work yet    
#@skip    
def test_director_serialization(isolated_dir):
    '''Test where the ILS is a python derived from director class'''
    ils = IlsInPython()
    res = ils.apply_ils([1,2,3], [4,5,6], [0,1,2])
    assert np.allclose(res, np.array([1.0,2.0,3.0]))
    write_shelve("temp.xml", ils)
    # Note serialize_write_string doesn't work. The pickle object contains
    # items that aren't utf-8. May be some way to work around this (e.g.,
    # return bytes instead of string), but it isn't clear that this is all
    # that important. We can save to a file, and we can save to binary
    ilsr = read_shelve("temp.xml")
    res = ilsr.apply_ils([1,2,3], [4,5,6], [0,1,2])
    assert np.allclose(res, np.array([1.0,2.0,3.0]))
