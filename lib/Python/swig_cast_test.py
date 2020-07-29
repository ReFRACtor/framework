from builtins import object
# Test that we do casting correctly, where we take a function returning a
# boost::shared_ptr<Base> and return the exact derived type to python
# (e.g., a function returning a boost::shared_ptr<Pressure>, which is actually
# a PressureSigma, will return the PressureSigma type to python.

from .test_support import *

# A simple python based class, to 
class PythonPressureSigma(rf.PressureImpBase):
    '''This class uses pressure sigma levels to determine the pressure
    levels we do the RT on.'''
    def __init__(self, a, b, surface_pressure, pressure_flag = True):
        coef = np.array([surface_pressure])
        flag = np.array([pressure_flag])
        rf.PressureImpBase.__init__(self,coef, flag)
        self.a = a
        self.b = b

    def clone(self):
        res = rf.PressureSigma(self.a, self.b,
                            self.coefficient().value()[0],
                            self.used_flag_value()[0])
        return res
    
    def calc_pressure_grid(self):
        t = self.b * self.coefficient()[0] + self.a
        self.pgrid.reference(np_to_array_ad(t))

    def state_vector_name_i(self, i):
        return "Surface Pressure (Pascals)"

    def my_func(self):
        '''Function just to make sure we have an object this type.'''
        return 102

    def desc(self):
        return '''
Pressure Sigma:
  Surface pressure: %f
  Retrieval flag: %s
  a: %s
  b: %s
''' %(self.coefficient().value()[0], self.used_flag_value()[0].__str__(),
      self.a.__str__(), self.b.__str__())

def test_cast_cpp():
    '''Test returning a C++ class'''
    psigma = rf.PressureSigma([0,0,0], [0.3, 0.6, 1.0], 10, True)
    pwrap = rf.PressureHolder(psigma)
    # Test functions only in PressureSigma
    assert pwrap.p.b == approx([0.3, 0.6, 1.0])

    pwrap.p = rf.PressureSigma([0,0,0], [0.4, 0.7, 1.1], 10, True)
    assert pwrap.p.b == approx([0.4, 0.7, 1.1])

def test_cast_cpp_excep():
    '''Test returning a C++ class'''
    psigma = rf.PressureSigma([0,0,0], [0.3, 0.6, 1.0], 10, True)
    pwrap = rf.PressureHolder(psigma)
    with pytest.raises(AttributeError):
        pwrap.p.temperature

def test_cast_python():
    '''Make sure we handle classes that are actually python correctly.'''
    psigma = PythonPressureSigma([0,0,0], [0.3, 0.6, 1.0], 10, True)
    pwrap = rf.PressureHolder(psigma)
    assert pwrap.p.my_func() == 102
