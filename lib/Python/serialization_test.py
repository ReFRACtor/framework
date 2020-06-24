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

# We have a few things that we need to test in a lot of places. Not sure
# we want this equality everywhere, this is really just "true for unit
# tests." We may move these into a more central place like test_support.py,
# but for now have here.

def _double_with_unit_eq(self, v):
    return (self.value == pytest.approx(v.value)
            and self.units == v.units)

DoubleWithUnit.__eq__ = _double_with_unit_eq

def _spectral_domain_eq(s1, s2):
    # Let assertion exceptions pass through. This is kind of breaking
    # what == is suppose to be, but since these end up in a assert statement
    # it is useful to trace back to what didn't match.
    assert_allclose(s1.data, s2.data)
    assert_allclose(s1.sample_index, s2.sample_index)
    assert s1.units == s2.units
    return True

SpectralDomain.__eq__ = _spectral_domain_eq

def _spectral_range_eq(s1, s2):
    # Let assertion exceptions pass through. This is kind of breaking
    # what == is suppose to be, but since these end up in a assert statement
    # it is useful to trace back to what didn't match.
    assert_allclose(s1.data, s2.data)
    assert_allclose(s1.uncertainty, s2.uncertainty)
    assert s1.units == s2.units
    return True

SpectralRange.__eq__ = _spectral_range_eq

def _array_with_unit_eq(a1, a2):
    assert_allclose(a1.value, a2.value)
    assert a1.units == a2.units
    return True

def _array_ad_eq(a1, a2):
    assert_allclose(a1.value, a2.value)
    assert_allclose(a1.jacobian, a2.jacobian)
    return True

def _array_ad_with_unit_eq(a1, a2):
    a1.value == a2.value
    assert a1.units == a2.units
    return True

ArrayWithUnit_double_1.__eq__ = _array_with_unit_eq
ArrayWithUnit_double_2.__eq__ = _array_with_unit_eq
ArrayWithUnit_double_3.__eq__ = _array_with_unit_eq
ArrayAd_double_1.__eq__ = _array_ad_eq
ArrayAd_double_2.__eq__ = _array_ad_eq
ArrayAd_double_3.__eq__ = _array_ad_eq
ArrayAdWithUnit_double_1.__eq__ = _array_ad_with_unit_eq
ArrayAdWithUnit_double_2.__eq__ = _array_ad_with_unit_eq
ArrayAdWithUnit_double_3.__eq__ = _array_ad_with_unit_eq

# Add auto_derivative

class TestPressureSigma(BaseForTesting):
    def create(self):
        # This is what is used in C++ unit tests in/meteorology/pressure_in
        pdata = [ 1.0000000000e+02, 7.0000000000e+03, 1.0000000000e+04,
          2.0000000000e+04, 2.8000000000e+04, 3.5000000000e+04,
          4.0000000000e+04, 4.5000000000e+04, 5.0000000000e+04,
          5.5000000000e+04, 6.0000000000e+04, 6.5000000000e+04,
          7.0000000000e+04, 7.5000000000e+04, 8.0000000000e+04,
          8.5000000000e+04, 9.0000000000e+04, 9.5000000000e+04,
          9.6716624900e+04  ]
        return PressureSigma(pdata, pdata[-1], False)

    def check(self, psigma, psigma2):
        assert_allclose(psigma.a, psigma2.a)
        assert_allclose(psigma.b, psigma2.b)
        
class TestArrayAd(BaseForTesting):
    def create(self):
        return ArrayAd_double_1([1,2,3],[[1,0,0],[0,1,0],[0,0,1]])

    def check(self, a1, a2):
        assert a1 == a2

class TestArrayAdWithUnit(BaseForTesting):
    def create(self):
        self.t = TestArrayAd()
        return ArrayAdWithUnit_double_1(self.t.create(), Unit("m/s"))

    def check(self, a1, a2):
        assert a1 == a2

class TestArrayWithUnit(BaseForTesting):
    def create(self):
        return ArrayWithUnit_double_1([1, 2, 3], Unit("m/s"))

    def check(self, a1, a2):
        a1 == a2

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
        assert v1 == v2
        
class TestDefaultConstant(BaseForTesting):
    def create(self):
        return DefaultConstant()

    def check(self, c1, c2):
        assert (c1.rayleigh_depolarization_factor ==
                pytest.approx(c2.rayleigh_depolarization_factor))
        assert c1.rayleigh_a == c2.rayleigh_a
        assert c1.rayleigh_b == c2.rayleigh_b
        assert c1.molar_weight_dry_air == c2.molar_weight_dry_air
        assert c1.molar_weight_water == c2.molar_weight_water
        assert c1.avogadro_constant == c2.avogadro_constant

class TestHdfFile(BaseForTesting):
    def create(self):
        return HdfFile(unit_test_data + "/in/common/l1b_example_data.h5")

    def check(self, h1, h2):
        assert h1.file_name == h2.file_name
        assert h1.mode == h2.mode
        assert_allclose(h1.read_double_3d("/Level1b/stokes_coefficient"),
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
        for i in range(b1.number_spectrometer):
            assert b1.lower_bound(i) == b2.lower_bound(i)
            assert b1.upper_bound(i) == b2.upper_bound(i)

class TestSpectralDomain(BaseForTesting):
    def create(self):
        return SpectralDomain([1,2,3], Unit("micron"))

    def check(self, s1, s2):
        assert s1 == s2
        
class TestSpectralRange(BaseForTesting):
    def create(self):
        return SpectralRange([1,2,3], Unit("m/s"), [0.1,0.2,0.3])

    def check(self, s1, s2):
        assert s1 == s2

class TestSpectrum(BaseForTesting):
    def create(self):
        self.t1 = TestSpectralDomain()
        self.t2 = TestSpectralRange()
        return Spectrum(self.t1.create(), self.t2.create())
    
    def check(self, s1, s2):
        assert s1.spectral_domain == s2.spectral_domain
        assert s1.spectral_range == s2.spectral_range

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
        assert_allclose(p1.pressure_level, p2.pressure_level)

class TestPressureFixedLevel(BaseForTesting):
    def create(self):
        self.t = TestPressureLevelInput()
        return PressureFixedLevel(True, self.t.create(), 3)
    
    def check(self, p1, p2):
        assert p1.pressure_grid == p2.pressure_grid

class TestTemperatureFixedLevel(BaseForTesting):
    def create(self):
        self.t = TestPressureLevelInput()
        self.t2 = TestPressureFixedLevel()
        return TemperatureFixedLevel([True,True,False],False,[10,11,12],
                                     0, self.t2.create(), self.t.create())

    def check(self, temp1, temp2):
        pres = self.t2.create()
        assert temp1.temperature_grid(pres) == temp2.temperature_grid(pres)

class TestAbsorberVmrFixedLevel(BaseForTesting):
    def create(self):
        self.t = TestPressureLevelInput()
        self.t2 = TestPressureFixedLevel()
        return AbsorberVmrFixedLevel(self.t2.create(),
                                     self.t.create(),
                                     [True,True,False],[10,11,12],"Fred")

    def check(self, a1, a2):
        assert a1.gas_name == a2.gas_name
        assert_allclose(a1.volume_mixing_ratio_level,
                        a2.volume_mixing_ratio_level)
        press = self.t2.create()
        assert a1.vmr_grid(press) == a2.vmr_grid(press)

class TestAbsorberVmrFixedLevelScaled(BaseForTesting):
    def create(self):
        self.t = TestPressureLevelInput()
        self.t2 = TestPressureFixedLevel()
        return AbsorberVmrFixedLevelScaled(self.t2.create(),
                                           self.t.create(),
                                           [10,11,12], True,
                                           1.5, "Fred")

    def check(self, a1, a2):
        assert a1.gas_name == a2.gas_name
        assert a1.scale_factor == pytest.approx(a2.scale_factor)
        press = self.t2.create()
        assert a1.vmr_grid(press) == a2.vmr_grid(press)

class TestExampleLevelL1b(BaseForTesting):
    def create(self):
        self.t = TestHdfFile()
        return ExampleLevel1b(self.t.create(), "2014090915251774")

    def check(self, f1, f2):
        self.t.check(f1.input, f2.input)
        assert f1.data_index == f2.data_index

class TestLevelL1bCache(BaseForTesting):
    def create(self):
        self.t = TestHdfFile()
        return Level1bCache(ExampleLevel1b(self.t.create(), "2014090915251774"))

    def check(self, f1, f2):
        assert f1.number_spectrometer == f2.number_spectrometer
        for i in range(f1.number_spectrometer):
            assert f1.latitude(i) == f2.latitude(i)
            assert f1.longitude(i) == f2.longitude(i)
            assert f1.sounding_zenith(i) == f2.sounding_zenith(i)
            assert f1.sounding_azimuth(i) == f2.sounding_azimuth(i)
            assert f1.solar_zenith(i) == f2.solar_zenith(i)
            assert f1.solar_azimuth(i) == f2.solar_azimuth(i)
            assert f1.altitude(i) == f2.altitude(i)
            assert f1.relative_velocity(i) == f2.relative_velocity(i)
            assert f1.sample_grid(i) == f2.sample_grid(i)
            assert f1.radiance(i) == f2.radiance(i)
            assert f1.time(i).pgs_time == pytest.approx(f2.time(i).pgs_time)
        
class TestExampleLevelL1bInfo(BaseForTesting):
    def create(self):
        self.t = TestHdfFile()
        return ExampleLevel1bInfo(self.t.create())

    def check(self, f1, f2):
        self.t.check(f1.input, f2.input)
        lst1 = f1.level1b_list()
        lst2 = f2.level1b_list()
        t2 = TestExampleLevelL1b()
        t2.create()
        assert len(lst1) == len(lst2)
        for i in range(len(lst1)):
            t2.check(GenericObject.convert_to_most_specific_class(lst1[i]),
                     GenericObject.convert_to_most_specific_class(lst2[i]))
        

class TestExampleMetFile(BaseForTesting):
    def create(self):
        return ExampleMetFile(unit_test_data +
                              "in/meteorology/example_met_data.h5",
                              "20091009203401")

    def check(self, f1, f2):
        assert f1.input.file_name == f2.input.file_name
        assert f1.input.mode == f2.input.mode
        assert f1.data_index == f2.data_index
        assert_allclose(f1.specific_humidity, f2.specific_humidity)
            
class TestTemperatureLevelOffset(BaseForTesting):
    def create(self):
        self.t = TestPressureSigma()
        temp_d = [244.2, 214.553, 218.029, 222.544, 218.341, 221.37,
                  227.38, 233.493, 239.376, 244.52, 248.708, 251.979,
                  254.537, 256.655, 258.521, 260.155, 261.747,
                  261.732, 258.598]
        return TemperatureLevelOffset(self.t.create(), temp_d, 0, True)

    def check(self, temp1, temp2):
        pres = self.t.create()
        assert temp1.temperature_grid(pres) == temp2.temperature_grid(pres)

class TestTemperatureMet(BaseForTesting):
    def create(self):
        self.t = TestPressureSigma()
        self.t2 = TestExampleMetFile()
        return TemperatureMet(self.t2.create(), self.t.create(), 0, True)

    def check(self, temp1, temp2):
        pres = self.t.create()
        assert temp1.temperature_grid(pres) == temp2.temperature_grid(pres)

class TestAltitudeHydrostatic(BaseForTesting):
    def create(self):
        self.t = TestPressureSigma()
        self.t2 = TestTemperatureMet()
        return AltitudeHydrostatic(self.t.create(), self.t2.create(),
                                   DoubleWithUnit(77.1828918457, "deg"),
                                   DoubleWithUnit(416, "m"))

    def check(self, a1, a2):
        pgrid = self.t.create().pressure_grid
        t3 = TestAutoDerivativeWithUnit()
        t3.create()
        for i in range(pgrid.value.rows):
            p_i = AutoDerivativeWithUnitDouble(pgrid.value[i], pgrid.units)
            t3.check(a1.altitude(p_i), a2.altitude(p_i))
            t3.check(a1.gravity(p_i), a2.gravity(p_i))
        
class TestSurfaceTemperatureDirect(BaseForTesting):
    def create(self):
        return SurfaceTemperatureDirect(ArrayWithUnit_double_1([258.598,259.598], "K"),
                                        [False, False])

    def check(self, t1, t2):
        t = TestAutoDerivativeWithUnit()
        t.create()
        t.check(t1.surface_temperature(0), t2.surface_temperature(0))
        t.check(t1.surface_temperature(1), t2.surface_temperature(1))
        
class TestRayleighBodhaine(BaseForTesting):
    def create(self):
        self.t = TestDefaultConstant()
        self.t2 = TestPressureSigma()
        self.t3 = TestAltitudeHydrostatic()
        return RayleighBodhaine(self.t2.create(), [self.t3.create(),],
                                self.t.create())

    def check(self, r1, r2):
        wl = DoubleWithUnit(3.091472458622220074e+02, "nm")
        r1.cross_section(wl) == r2.cross_section(wl)

class TestRayleighYoung(BaseForTesting):
    def create(self):
        self.t = TestDefaultConstant()
        self.t2 = TestPressureSigma()
        self.t3 = TestAltitudeHydrostatic()
        return RayleighYoung(self.t2.create(), [self.t3.create(),],
                                self.t.create())

    def check(self, r1, r2):
        wl = DoubleWithUnit(3.091472458622220074e+02, "nm")
        r1.cross_section(wl) == r2.cross_section(wl)

class TestCompositeInitialGuess(BaseForTesting):
    def create(self):
        ig = InitialGuessValue()
        ig.apriori = [1,2]
        ig.apriori_covariance = [[1,2],[3,4]]
        ig2 = InitialGuessValue()
        ig.apriori = [1,2, 3]
        ig.apriori_covariance = [[1,2, 3],[4,5,6],[7,8,9]]
        ci = CompositeInitialGuess()
        ci.add_builder(ig)
        ci.add_builder(ig2)
        return ci
    
    def check(self, ig1, ig2):
        assert_allclose(ig1.initial_guess, ig2.initial_guess)
        assert_allclose(ig1.apriori, ig2.apriori)
        assert_allclose(ig1.apriori_covariance, ig2.apriori_covariance)

class TestBardNLSSProblem(BaseForTesting):
    def create(self):
        f = BardNLLSProblem()
        f.parameters = [1.0,1.0,1.0]
        return f

    def check(self, f1, f2):
        assert_allclose(f1.residual, f2.residual)

class TestMeyerNLSSProblem(BaseForTesting):
    def create(self):
        f = MeyerNLLSProblem()
        f.parameters = [0.02, 4000.0, 250.0]
        return f

    def check(self, f1, f2):
        assert_allclose(f1.residual, f2.residual)

class TestBardMLProblem(BaseForTesting):
    def create(self):
        # From the nlls_solver_gsl unit test
        measurement = [0.14, 0.18, 0.22, 0.25, 0.29, 0.32,
                       0.35, 0.39, 0.37, 0.58, 0.73, 0.96, 1.34, 2.10, 4.39]
        measurement_error_cov = np.full_like(measurement, 1.0)
        
        ml = BardMLProblem(measurement, measurement_error_cov)
        ml.parameters = [1.0,1.0,1.0]
        return ml

    def check(self, ml1, ml2):
        assert_allclose(ml1.model, ml2.model)

class TestMeyerMLProblem(BaseForTesting):
    def create(self):
        # From the nlls_solver_gsl unit test
        measurement = [34780.0, 28610.0, 23650.0, 19630.0, 16370.0,
                       13720.0, 11540.0, 9744.0, 8261.0, 7030.0, 6005.0,
                       5147.0, 4427.0, 3820.0, 3307.0, 2872.0]
        measurement_error_cov = np.full_like(measurement, 1.0)
        
        ml = MeyerMLProblem(measurement, measurement_error_cov)
        ml.parameters = [0.02, 4000.0, 250.0]
        return ml

    def check(self, ml1, ml2):
        assert_allclose(ml1.model, ml2.model)
        
class TestCostMinimizerGSL(BaseForTesting):
    def create(self):
        # We save 2 of these, both at the start and after solving.
        # We then use the start version from serialization, run, and
        # make sure it matches the finished version.
        # The data comes from the cost_minimizer_gsl_test.cc unit test
        m = GenericObjectMap()
        x0 = [1.0,1.0,1.0]
        p1 = BardNLLSProblem()
        p1.parameters = x0
        solv1 = CostMinimizerGSL(p1, 1000, 1e-5)
        p2 = BardNLLSProblem()
        p2.parameters = x0
        solv2 = CostMinimizerGSL(p2, 1000, 1e-5)
        solv2.solve()
        m["solv_init"] = solv1
        m["solv_after"] = solv2
        return m

    def compare(self, s1, s2):
        '''Compare two solvers'''
        assert s1.status == s2.status
        assert s1.num_accepted_steps == s2.num_accepted_steps
        ap1 = s1.accepted_points
        ap2 = s2.accepted_points
        for i in range(len(ap1)):
            assert_almost_equal(ap1[i], ap2[i])
        assert_almost_equal(s1.cost_at_accepted_points,
                            s2.cost_at_accepted_points)
        
    def check(self, m1, m2):
        self.compare(m1.solv_init, m2.solv_init)
        self.compare(m1.solv_after, m2.solv_after)
        solv = m2.solv_init
        solv.solve()
        self.compare(m1.solv_after, solv)

class BaseForIterativeSolverDer(BaseForTesting):
    '''There a a few solvers that derive from IterativeSolverDer. These
    can all be tested the same way, although the creation will be 
    different. Pull out the comparison parts'''
    def create_problem(self):
        x0 = [1.0,1.0,1.0]
        p = BardNLLSProblem()
        p.parameters = x0
        return p
        
    def compare(self, s1, s2):
        '''Compare two solvers'''
        assert s1.status == s2.status
        assert s1.num_accepted_steps == s2.num_accepted_steps
        ap1 = s1.accepted_points
        ap2 = s2.accepted_points
        for i in range(len(ap1)):
            assert_almost_equal(ap1[i], ap2[i])
        gap1 = s1.gradient_at_accepted_points
        gap2 = s2.gradient_at_accepted_points
        for i in range(len(gap1)):
            assert_almost_equal(gap1[i], gap2[i])
        assert_almost_equal(s1.cost_at_accepted_points,
                            s2.cost_at_accepted_points)
        
    def check(self, m1, m2):
        self.compare(m1.solv_init, m2.solv_init)
        self.compare(m1.solv_after, m2.solv_after)
        solv = m2.solv_init
        solv.solve()
        self.compare(m1.solv_after, solv)

class TestNLLSSolverGSLSM(BaseForIterativeSolverDer):
    def create(self):
        m = GenericObjectMap()
        solv1 = NLLSSolverGSLSM(self.create_problem(), 100)
        solv2 = NLLSSolverGSLSM(self.create_problem(), 100)
        solv2.solve()
        m["solv_init"] = solv1
        m["solv_after"] = solv2
        return m

class TestNLLSSolverGSL(BaseForIterativeSolverDer):
    def create(self):
        m = GenericObjectMap()
        solv1 = NLLSSolverGSL(self.create_problem(), 100)
        solv2 = NLLSSolverGSL(self.create_problem(), 100)
        solv2.solve()
        m["solv_init"] = solv1
        m["solv_after"] = solv2
        return m

class TestRelativeHumidity(object):
    def test_serialization(self, sample_absorber, sample_temperature,
                           sample_pressure):
        v1 = RelativeHumidity(sample_absorber, sample_temperature,
                              sample_pressure);
        t = serialize_write_string(v1)
        logging.info("Serialization:\n%s" % t)
        v2 = serialize_read_generic_string(t)
        self.check(v1, v2)

    def check(self, v1, v2):
        return
        assert v1.specific_humidity_grid == v2.specific_humidity_grid
        assert v1.relative_humidity_grid == v2.relative_humidity_grid
        assert v1.relative_humidity_layer == v2.relative_humidity_layer

class TestAerosolExtinctionLinear(BaseForTesting):
    def create(self):
        self.t = TestPressureSigma()
        return AerosolExtinctionLinear(self.t.create(), [True, True, True],
                                       [0.0, 1.0, 0.2], "Kahn")

    def check(self, a1, a2):
        assert a1.aerosol_extinction == a2.aerosol_extinction
        assert a1.aerosol_name == a2.aerosol_name
        assert a1.model_short_name == a2.model_short_name
