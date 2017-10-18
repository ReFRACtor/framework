import os
import h5py
import logging

import numpy as np

import refractor.factory.creator as creator
import refractor.factory.param as param
from refractor.factory import process_config
from refractor import framework as rf

logging.basicConfig(level=logging.DEBUG)

static_input_file = os.path.join(os.path.dirname(__file__), "../lua/example_static_input.h5")
static_input = h5py.File(static_input_file)

data_dir = os.path.join(os.path.dirname(__file__), '../in/common')
l1b_file = os.path.join(data_dir, "l1b_example_data.h5")
met_file = os.path.join(data_dir, "met_example_data.h5")

observation_id = "2014090915251774"

# Helpers to abstract away getting data out of the static input file
def static_value(dataset):
    return static_input[dataset][:]

def static_units(dataset):
    return static_input[dataset].attrs['Units'][0].decode('UTF8') 


config_def = {
    'creator': creator.base.SaveToCommon,
    'order': ['common', 'input', 'spec_win', 'spectrum_sampling', 'atmosphere'],
    'common': {
        'creator': creator.base.SaveToCommon,
        'desc_band_name': static_value("Common/desc_band_name"),
        'hdf_band_name': static_value("Common/hdf_band_name"),
        'band_reference': {
            'creator': creator.value.ArrayWithUnit,
            'value': static_value("Common/band_reference_point"),
            'units': static_units("Common/band_reference_point"),
        },
        'num_channels': 3,
        'absco_base_path': '/mnt/data1/absco',
        'constants': {
            'creator': creator.common.DefaultConstants,
        },
    },
    'input': {
        'creator': creator.base.SaveToCommon,
        'l1b': rf.ExampleLevel1b(l1b_file, observation_id),
        'met': rf.ExampleMetFile(met_file, observation_id),
    },
    'spec_win': {
        'creator': creator.forward_model.SpectralWindowRange,
        'window_ranges': {
            'creator': creator.value.ArrayWithUnit,
            'value': static_value("/Spectral_Window/microwindow"),
            'units': static_units("/Spectral_Window/microwindow"),
        },
    },
    'spectrum_sampling': {
        'creator': creator.forward_model.UniformSpectrumSampling,
        'high_res_spacing': rf.DoubleWithUnit(0.01, "cm^-1"), 
    },
    'stokes_coefficient': {
    },
    'instrument': {
    },
    'spectrum_effect': {
    },
    'rt': {
    },
    'state_vector': {
    },
    'atmosphere': {
        'creator': creator.atmosphere.AtmosphereCreator,
        'pressure': {
            'creator': creator.atmosphere.PressureSigma,
            'apriori': {
                'creator': creator.met.ValueFromMet,
                'field': "surface_pressure",
            },
            'a_coeff': static_value("Pressure/Pressure_sigma_a"),
            'b_coeff': static_value("Pressure/Pressure_sigma_b"),
        },
        'temperature': {
            'creator': creator.atmosphere.TemperatureMet,
            'apriori': static_value("Temperature/Offset/a_priori")
        },
        'altitudes': { 
            'creator': creator.atmosphere.AltitudeHydrostatic,
            'latitude': {
                'creator': creator.l1b.ValueFromLevel1b,
                'field': "latitude",
            },
            'surface_height': {
                'creator': creator.l1b.ValueFromLevel1b,
                'field': "altitude",
            },
        },
        'absorber': {
            'creator': creator.atmosphere.AbsorberAbsco,
            'gases': ['CO2'],
            'CO2': {
                'creator': creator.atmosphere.AbsorberGasDefinition,
                'vmr': {
                    'creator': creator.atmosphere.AbsorberVmrLevel,
                    'apriori': {
                        'creator': creator.atmosphere.GasVmrApriori,
                        'gas_name': 'CO2',
                        'reference_atm_file': static_input_file,
                    },
                },
                'absorption': {
                    'creator': creator.atmosphere.AbscoHdf,
                    'table_scale': [1.0, 1.0, 1.004],
                    'filename': "v5.0.0/co2_devi2015_wco2scale-nist_sco2scale-unity.h5",
                },
            },
            'H2O': {
            },
            'O2': {
            },
        },
    },
}

config_inst = process_config(config_def)

#print(config_inst
from pprint import pprint
pprint(config_inst, indent=4)

