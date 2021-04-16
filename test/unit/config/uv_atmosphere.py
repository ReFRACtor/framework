import os
import re
from copy import deepcopy

import numpy as np

import refractor.factory.creator as creator
from refractor import framework as rf
from refractor.config import refractor_config

cross_section_dir = os.path.join(os.path.dirname(__file__), '../../../input/cross_sections')
profile_filename = os.path.join(os.path.dirname(__file__), '../data/in/uv_atmosphere/Profiles_9_2006726_1500.dat')

cross_section_filenames = {
    'O3': "o3abs_brion_195_660_vacfinal.dat", 
    'NO2': "no2r_97.nm",
    'SO2': "SO2_298_BISA.nm",
    'HCHO': "H2CO_Meller_Moortgat_MPI.txt",
    'BrO': "228kbro10cm_padded.nm",
    'CHOCHO': "CHOCHO_Xsections_250-510nm.dat",
}

# Add path to filenames
for gas_name, filename in cross_section_filenames.items():
    cross_section_filenames[gas_name] = os.path.join(cross_section_dir, filename)
 
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

def vmr_profile(gas_name):
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
                        'vmr_profile': None,
                    },
                    'cross_section': {
                        'creator': creator.absorber.CrossSectionTableAscii,
                        'filename': None,
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

    # Create cross section table set ups for each gass
    for gas_name in config_def['atmosphere']['absorber']['gases']:
        gas_def = deepcopy(config_def['atmosphere']['absorber']['default_gas_definition'])

        gas_def['vmr']['vmr_profile'] = vmr_profile(gas_name)
        gas_def['cross_section']['filename'] = cross_section_filenames[gas_name]

        config_def['atmosphere']['absorber'][gas_name] = gas_def

    # Set the O3 conversion factor
    config_def['atmosphere']['absorber']['O3']['cross_section']['conversion_factor'] = 1e20
    config_def['atmosphere']['absorber']['CHOCHO']['cross_section']['conversion_factor'] = 5e18

    return config_def
