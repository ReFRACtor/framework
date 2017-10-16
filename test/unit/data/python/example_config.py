import os
import h5py
import logging

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

config_def = {
    'order': ['common', 'input', 'spec_win', 'spectrum_sampling'],
    'common': {
        'creator': creator.base.SaveToCommon,
        'desc_band_name': static_input["Common/desc_band_name"][:],
        'hdf_band_name': static_input["Common/hdf_band_name"][:],
        'band_reference': {
            'creator': creator.value.ArrayWithUnit,
            'value': static_input["Common/band_reference_point"][:],
            'units': static_input["Common/band_reference_point"].attrs['Units'][0].decode('UTF8'),
        },
        'num_channels': 4,
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
            'value': static_input["/Spectral_Window/microwindow"][:],
            'units': static_input["/Spectral_Window/microwindow"].attrs['Units'][0].decode('UTF8')
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
    },
}

config_inst = process_config(config_def)

#print(config_inst
from pprint import pprint
pprint(config_inst, indent=4)

