import os
import re
from copy import deepcopy

import numpy as np

from refractor import framework as rf
from refractor.framework import creator, refractor_config
from refractor.framework import cross_section_filenames, cross_section_file_conversion_factors

profile_filename = os.path.join(os.path.dirname(__file__), '../data/in/uv_atmosphere/Profiles_9_2006726_1500.dat')
 
class ProfileReader(object):

    def __init__(self, filename):
        self.header = {}

        with open(filename) as prof_file:
            # Skip first line
            line = ""
            while line.find("End_of_Header") < 0:
                line = prof_file.readline()
                if line.find("=") >= 0:
                    k, v = re.split('\s+=\s+', line.strip())
                    self.header[k] = v

            self.column_names = prof_file.readline().strip().split()
            self.unit_names = prof_file.readline().strip().split()

            self.data = np.loadtxt(prof_file)

        self.latitude = rf.ArrayWithUnit(np.array([float(self.header['latitude'])]), "deg")
        self.longitude = rf.ArrayWithUnit(np.array([float(self.header['longitude'])]), "deg")
        self.surface_height = rf.ArrayWithUnit(np.array([float(self.header['ZSUR(m)'])]), "m")

profile = ProfileReader(profile_filename)

def surface_pressure():
    press_col = profile.column_names.index("Pressure")
    # hPa -> Pa
    return profile.data[0, press_col] * 100

def pressure_levels():
    press_col = profile.column_names.index("Pressure")
    # hPa -> Pa
    return profile.data[::-1, press_col] * 100

def temperature_levels():
    temp_col = profile.column_names.index("TATM")
    return profile.data[::-1, temp_col]

def vmr_profile(gas_name=None):
    if gas_name == "CHOCHO":
        col_name = "CHOHO"
    else:
        col_name = gas_name
    vmr_col = profile.column_names.index(col_name)
    return profile.data[::-1, vmr_col]
       
# Common configuration definition shared amongst retrieval and simulation types of configuration
@refractor_config
def uv_atmosphere_definition():

    config_def = {
        'creator': creator.base.ParamPassThru,
        'order': ['common', 'scene', 'atmosphere'],
        'child': 'atmosphere',
        'common': {
            'creator': creator.base.SaveToCommon,
            'band_reference': rf.ArrayWithUnit(np.array([300.0]), "nm"),
            'desc_band_name': ['UV'],
            'constants': {
                'creator': creator.common.DefaultConstants,
            },
        },
        'scene': {
            'creator': creator.base.SaveToCommon,
            'latitude': profile.latitude,
            'surface_height': profile.surface_height,
            'num_channels': 1,
        },
        'atmosphere': {
            'creator': creator.atmosphere.AtmosphereCreator,
            'pressure': {
                'creator': creator.atmosphere.PressureGrid,
                'surface_pressure': surface_pressure,
                'pressure_levels': pressure_levels,
            },
            'temperature': {
                'creator': creator.atmosphere.TemperatureLevel,
                'temperature_profile': temperature_levels,
            },
            'altitude': {
                'creator': creator.atmosphere.AltitudeHydrostatic,
            },
            'rayleigh': {
                'creator': creator.rayleigh.RayleighBodhaine,
            },
            'absorber': {
                'creator': creator.absorber.AbsorberXSec,
                'gases': ['O3', 'NO2', 'SO2', 'HCHO', 'BrO', 'CHOCHO'],
                'default_gas_definition': {
                    'creator': creator.absorber.CrossSectionGasDefinition,
                    'vmr': {
                        'creator': creator.absorber.AbsorberVmrLevel,
                        'vmr_profile': vmr_profile,
                    },
                    'cross_section': {
                        'creator': creator.absorber.CrossSectionTableAscii,
                        'filename': None,
                        'filename': lambda gas_name=None: cross_section_filenames[gas_name],
                        'conversion_factor': lambda gas_name=None: cross_section_file_conversion_factors.get(gas_name, 1.0),
                    },
                }, 
            },
            'relative_humidity': {
                'creator': creator.atmosphere.RelativeHumidity,
            },
            'ground': {
                'creator': creator.base.PickChild,
                'child': 'lambertian',
                'lambertian': {
                    'creator': creator.ground.GroundLambertian,
                    'polynomial_coeffs': np.array([[0.1, 0]]),
                },
            },
        },
    }

    return config_def
