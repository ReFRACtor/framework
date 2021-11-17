from refractor.framework import *
import pickle
import h5py

fmodel_from_lua = pickle.load(open("fmodel.pkl", "rb"))
atm_from_lua = fmodel_from_lua.radiative_transfer.low_stream_radiative_transfer.atmosphere

sid = "2014090915251774"
base_dir = "/home/smyth/Local/framework/"

class my_cached_property:
    '''Python 3.8 introduced a cached_property decorator in functools.
    That is a little new to depend on, so we'll write our own simple
    one.

    Note that caching here is not just about the slight performance
    advantage.  The real reason we cache is that we need the same
    object returned every time. So if we grab a Pressure to use in an
    AerosolProperty, we need the exact same Pressure object also used
    to give to a Temperature.  This is important because the Pressure
    might get updated by a StateVector change, and we need
    AerosolProperty and Temperature to see the same update.

    This class can go away with python 3.8 is old enough that we can depend
    on using it and having the system cached_property available.
    '''
    def __init__(self, func):
        self.func = func
        self.attr_name = "_" + func.__name__

    def __get__(self, instance, objtype):
        if(not hasattr(instance, self.attr_name)):
            setattr(instance, self.attr_name, self.func(instance))
        return getattr(instance, self.attr_name)

class ConfigSolar:
    '''This has a few pieces, so we group it together'''
    def __init__(self, solar_file, l1b_file, spec_index):
        self.solar_file = solar_file
        self.l1b_file = l1b_file
        self.spec_index = spec_index
        
    @my_cached_property
    def solar_dopler_shift(self):
        return SolarDopplerShiftPolynomial \
            (self.l1b_file.time(self.spec_index),
             self.l1b_file.latitude(self.spec_index),
             self.l1b_file.solar_zenith(self.spec_index), 
             self.l1b_file.solar_azimuth(self.spec_index),
             self.l1b_file.altitude(self.spec_index),
             DefaultConstant(), True)

    @my_cached_property
    def solar_absorption(self):
        return SolarAbsorptionTable(self.solar_file,
                      "/Solar/Absorption/Absorption_%d" % (self.spec_index+1))
    
    @my_cached_property
    def solar_continuum(self):
        return SolarContinuumTable(self.solar_file,
                     "/Solar/Continuum/Continuum_%d" % (self.spec_index+1), False)

    @my_cached_property
    def solar_model(self):
        return SolarAbsorptionAndContinuum(self.solar_dopler_shift,
                                           self.solar_absorption,
                                           self.solar_continuum)
        
class BaseConfig:
    def __init__(self, sid, base_dir = "/home/smyth/Local/framework/"):
        self.static_file = h5py.File(base_dir + "/test/unit/data/lua/example_static_input.h5")
        self.met_file = ExampleMetFile(HdfFile(base_dir + "/test/unit/data/in/common/met_example_data.h5"), sid)
        self.l1b_file = ExampleLevel1b(HdfFile(base_dir + "/test/unit/data/in/common/l1b_example_data.h5"), sid)
        self.solar_file = HdfFile(base_dir + "input/common/input/l2_solar_model.h5")
        self.aerosol_file_name = base_dir + "input/common/input/l2_aerosol_combined.h5"
        self.desc_band_name = [t.decode('utf-8') for
                  t in self.static_file["Common/desc_band_name"][:]]
        self.hdf_band_name = [t.decode('utf-8') for
                  t in self.static_file["Common/hdf_band_name"][:]]
        self.band_reference = ArrayWithUnit_double_1(
            self.static_file["Common/band_reference_point"][:],
            self.static_file["Common/band_reference_point"].attrs['Units'][0].decode('utf-8'))
            
    @my_cached_property
    def pressure(self):
        a = self.static_file["Pressure/Pressure_sigma_a"][:]
        b = self.static_file["Pressure/Pressure_sigma_b"][:]
        return PressureSigma(a, b, self.met_file.surface_pressure, True)

    @my_cached_property
    def temperature(self):
        toff = self.static_file["Temperature/Offset/a_priori"][0]
        return TemperatureMet(self.met_file, self.pressure, toff, True)

    @my_cached_property
    def ground(self):
        # We'll come back to calculating ap. For now, just grab from Lua
        # result
        coeff = atm_from_lua.ground.coefficient.value
        ap = [[coeff[0],coeff[1]],
              [coeff[2],coeff[3]],
              [coeff[4],coeff[5]],]
        return GroundLambertian(ap, [[True,True],]*3,
                                self.band_reference, self.desc_band_name)

    @my_cached_property
    def spectral_window(self):
        win_range = ArrayWithUnit_double_3(
            self.static_file["Spectral_Window/microwindow"][:],
            self.static_file["Spectral_Window/microwindow"].attrs['Units'][0].decode('utf-8'))
        return SpectralWindowRange(win_range)

    @my_cached_property
    def solar_model(self):
        return [ConfigSolar(self.solar_file, self.l1b_file, i).solar_model
                for i in range(3)]
    
b = BaseConfig(sid)
print(b.pressure)
print(b.temperature)
print(b.ground)
print(b.spectral_window)
print(b.solar_model[0])
