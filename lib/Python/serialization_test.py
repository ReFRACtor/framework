from .test_support import *
from refractor_swig import *
import logging

# Note that we test some of the serialization at the C++ level with
# boost unit tests. However, a lot of this code is very cookie cutter,
# and it is just more convenient to test at python level. This also has
# the advantage that we test the python interface at the same time.
#
# We print out the serialization at the info level, you can view in
# pytest with --log-cli-level=info on the command line
#
# If you run into issues for a particular class, you may want to move the
# test down to the C++ level where you can cut out the overhead of building,
# installing, and running pytest. Debugging with pytest is certainly possible,
# but it is definitely easier at the C++ level. There are lots of examples
# of boost tests, see for example PressureSigma.
#
# Note the boost serialization is closely related to python pickling. For
# these tests we use our serialize_read_generic_string and
# serialize_write_string because we can somewhat view the xml generated.
# You can replace this with pickle dump and load - the only difference is you
# get a binary blob rather than xml test. For low level work, the xml can
# be read directly at the C++ level, while the pickle objects can't. But other
# than this, they pretty much operate the same.

class BaseForTesting(object):
    def create(self):
        pass

    def check(self, v1, v2):
        pass

    def test_serialization(self):
        v1 = self.create()
        t = serialize_write_string(v1)
        logging.info("Serialization:\n%s" % t)
        v2 = serialize_read_generic_string(t)
        self.check(v1, v2)

class TestPressureSigma(BaseForTesting):
    def create(self):
        return PressureSigma([0,0,0],[0.3,0.6,1.0], 10, True)

    def check(self, psigma, psigma2):
        assert_array_almost_equal(psigma.a, psigma2.a)
        assert_array_almost_equal(psigma.b, psigma2.b)
        
class TestArrayAd(BaseForTesting):
    def create(self):
        return ArrayAd_double_1([1,2,3],[[1,0,0],[0,1,0],[0,0,1]])

    def check(self, a1, a2):
        assert_array_almost_equal(a1.value, a2.value)
        assert_array_almost_equal(a1.jacobian, a2.jacobian)

class TestArrayAdWithUnit(BaseForTesting):
    def create(self):
        self.t = TestArrayAd()
        return ArrayAdWithUnit_double_1(self.t.create(), Unit("m/s"))

    def check(self, a1, a2):
        self.t.check(a1.value, a2.value)
        assert a1.units == a2.units

class TestArrayWithUnit(BaseForTesting):
    def create(self):
        return ArrayWithUnit_double_1([1, 2, 3], Unit("m/s"))

    def check(self, a1, a2):
        assert_array_almost_equal(a1.value, a2.value)
        assert a1.units == a2.units

class TestAutoDerivative(BaseForTesting):
    def create(self):
        return AutoDerivativeDouble(2, 0, 2)

    def check(self, v1, v2):
        assert v1.value == pytest.approx(v2.value)
        assert_array_almost_equal(v1.gradient, v2.gradient)

class TestAutoDerivativeWithUnit(BaseForTesting):
    def create(self):
        self.t = TestAutoDerivative()
        return AutoDerivativeWithUnitDouble(self.t.create(), Unit("m/s"))

    def check(self, a1, a2):
        self.t.check(a1.value, a2.value)
        assert a1.units == a2.units

class TestDoubleWithUnit(BaseForTesting):
    def create(self):
        return DoubleWithUnit(1.2345, Unit("m/s"))

    def check(self, v1, v2):
        assert v1.value == pytest.approx(v2.value)
        assert v1.units == v2.units
        
class TestDefaultConstant(BaseForTesting):
    def create(self):
        return DefaultConstant()

    def check(self, c1, c2):
        t = TestDoubleWithUnit()
        assert (c1.rayleigh_depolarization_factor ==
                pytest.approx(c2.rayleigh_depolarization_factor))
        t.check(c1.rayleigh_a, c2.rayleigh_a)
        t.check(c1.rayleigh_b, c2.rayleigh_b)
        t.check(c1.molar_weight_dry_air, c2.molar_weight_dry_air)
        t.check(c1.molar_weight_water, c2.molar_weight_water)
        t.check(c1.avogadro_constant, c2.avogadro_constant)

class TestHdfFile(BaseForTesting):
    def create(self):
        return HdfFile(unit_test_data + "/in/common/l1b_example_data.h5")

    def check(self, h1, h2):
        assert h1.file_name == h2.file_name
        assert h1.mode == h2.mode
        assert_array_almost_equal(
            h1.read_double_3d("/Level1b/stokes_coefficient"),
            h2.read_double_3d("/Level1b/stokes_coefficient"))

class TestMappingLinear(BaseForTesting):
    def create(self):
        return MappingLinear()

    def check(self, m1, m2):
        assert m1.name == m2.name

class TestMappingLog(BaseForTesting):
    def create(self):
        return MappingLog()

    def check(self, m1, m2):
        assert m1.name == m2.name

class TestMappingGaussian(BaseForTesting):
    def create(self):
        self.t = TestPressureSigma()
        return MappingGaussian(self.t.create(), True, 1e-10)

    def check(self, m1, m2):
        self.t.check(m1.pressure, m2.pressure)
        assert m1.name == m2.name
        assert m1.is_linear_total == m2.is_linear_total
        assert m1.min_desired == pytest.approx(m2.min_desired)

class TestMappingOffset(BaseForTesting):
    def create(self):
        return MappingOffset(1.0, [2, 3, 4])

    def check(self, m1, m2):
        assert m1.name == m2.name
        assert m1.initial_offset == pytest.approx(m2.initial_offset)
        assert_array_almost_equal(m1.offsetee, m2.offsetee)

class TestMappingScale(BaseForTesting):
    def create(self):
        return MappingScale(1.0, [2, 3, 4])

    def check(self, m1, m2):
        assert m1.name == m2.name
        assert m1.initial_scale_factor == pytest.approx(m2.initial_scale_factor)
        assert_array_almost_equal(m1.scalee, m2.scalee)

class TestSpectralBound(BaseForTesting):
    def create(self):
        return SpectralBound(ArrayWithUnit_double_2([[0.755, 0.785],
                                                     [1.58, 1.65],
                                                     [2.03, 2.09],
                                                     ],
                                                    Unit("micron")))

    def check(self, b1, b2):
        assert b1.number_spectrometer == b2.number_spectrometer
        t = TestDoubleWithUnit()
        for i in range(b1.number_spectrometer):
            t.check(b1.lower_bound(i), b2.lower_bound(i))
            t.check(b1.upper_bound(i), b2.upper_bound(i))

class TestSpectralDomain(BaseForTesting):
    def create(self):
        return SpectralDomain([1,2,3], Unit("micron"))

    def check(self, s1, s2):
        assert_array_almost_equal(s1.data, s2.data)
        assert_array_almost_equal(s1.sample_index, s2.sample_index)
        assert_array_almost_equal(s1.type_preference, s2.type_preference)
        assert s1.units == s2.units
        
class TestSpectralRange(BaseForTesting):
    def create(self):
        return SpectralRange([1,2,3], Unit("m/s"), [0.1,0.2,0.3])

    def check(self, s1, s2):
        assert_array_almost_equal(s1.data, s2.data)
        assert_array_almost_equal(s1.uncertainty, s2.uncertainty)
        assert s1.units == s2.units

class TestSpectrum(BaseForTesting):
    def create(self):
        self.t1 = TestSpectralDomain()
        self.t2 = TestSpectralRange()
        return Spectrum(self.t1.create(), self.t2.create())
    
    def check(self, s1, s2):
        self.t1.check(s1.spectral_domain, s2.spectral_domain)
        self.t2.check(s1.spectral_range, s2.spectral_range)

class TestNamedSpectrum(BaseForTesting):
    def create(self):
        self.t = TestSpectrum()
        return NamedSpectrum(self.t.create(), "Bob", 2);
    
    def check(self, s1, s2):
        self.t.check(s1, s2)
        assert s1.name == s2.name
        assert s1.index == s2.index

class TestStateVector(BaseForTesting):
    def create(self):
        sv = StateVector()
        sv.update_state([1,2,3], [[1,0.1,0.2],[0.1, 2, 0.3],[0.2,0.3,3]])
        return sv

    def check(self, s1, s2):
        assert_almost_equal(s1.state, s2.state)
        assert_almost_equal(s1.state_covariance, s2.state_covariance)

class TestPressureLevelInput(BaseForTesting):
    def create(self):
        return PressureLevelInput([1,2,3])

    def check(self, p1, p2):
        assert_almost_equal(p1.pressure_level, p2.pressure_level)

class TestPressureFixedLevel(BaseForTesting):
    def create(self):
        self.t = TestPressureLevelInput()
        return PressureFixedLevel(True, self.t.create(), 3)
    
    def check(self, p1, p2):
        t2 = TestArrayAdWithUnit()
        t2.create()
        t2.check(p1.pressure_grid, p2.pressure_grid)

class TestTemperatureFixedLevel(BaseForTesting):
    def create(self):
        self.t = TestPressureLevelInput()
        self.t2 = TestPressureFixedLevel()
        return TemperatureFixedLevel([True,True,False],False,[10,11,12],
                                     0, self.t2.create(), self.t.create())

    def check(self, temp1, temp2):
        t3 = TestArrayAdWithUnit()
        t3.create()
        pres = self.t2.create()
        t3.check(temp1.temperature_grid(pres), temp2.temperature_grid(pres))
    
        
