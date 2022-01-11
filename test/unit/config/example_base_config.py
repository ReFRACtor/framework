import os
import h5py
import logging

import numpy as np

from refractor.framework.factory import creator, param
from refractor.framework.config import refractor_config
from refractor import framework as rf

unit_test_dir = os.path.realpath(os.path.join(os.path.dirname(__file__), "../data"))
common_input_dir = os.path.realpath(os.path.join(os.path.dirname(__file__), "../../../input/common/input"))

static_input_file = os.path.join(unit_test_dir, "lua/example_static_input.h5")
ils_file = os.path.join(unit_test_dir, "lua/ils_data.h5")

solar_file = os.path.join(common_input_dir, "l2_solar_model.h5")
aerosol_prop_file = os.path.join(common_input_dir, "l2_aerosol_combined.h5")
covariance_file = os.path.join(os.path.dirname(__file__), "example_covariance.h5")

l1b_file = os.path.join(unit_test_dir, "in/common/l1b_example_data.h5")
met_file = os.path.join(unit_test_dir, "in/common/met_example_data.h5")

observation_id = "2014090915251774"
num_channels = 3

# Helpers to abstract away getting data out of the static input file
def static_value(dataset, dtype=None):
    with h5py.File(static_input_file, "r") as static_input:
        return np.array(static_input[dataset][:], dtype=dtype)

def static_scalar(dataset, dtype=None):
    value = static_value(dataset, dtype)
    if np.isscalar(value):
        return value
    else:
        if len(value) != 1:
            raise Exception(f"Expected scalar value, not array of size: {value.shape}")
        return value[0]

def static_units(dataset):
    with h5py.File(static_input_file, "r") as static_input:
        return static_input[dataset].attrs['Units'][0].decode('UTF8') 

def ils_delta_lambda():
    with h5py.File(ils_file) as ils_input:
        return ils_input["/InstrumentData/ils_delta_lambda"][:]

def ils_relative_response():
    with h5py.File(ils_file) as ils_input:
        return ils_input["/InstrumentData/ils_relative_response"][:]

class SampleGridIndexes(creator.Creator):

    instrument = param.InstanceOf(rf.Instrument)
    spec_win = param.InstanceOf(rf.SpectralWindow)
    spectrum_sampling = param.InstanceOf(rf.SpectrumSampling)
    
    num_channels = param.Scalar(int)

    def create(self, **kwargs):
        spectral_grid = rf.ForwardModelSpectralGrid(self.instrument(), self.spec_win(), self.spectrum_sampling())

        sample_indexes = []
        for chan_idx in range(self.num_channels()):
            chan_indexes = spectral_grid.pixel_list(chan_idx)
            
            # Remove last grid index so that unit tests are consistent with previous bug where
            # code was not being inclusive of the last point at the end of the supplied range
            chan_indexes = chan_indexes[:-1]

            sample_indexes.append(chan_indexes)

        return sample_indexes

@refractor_config
def base_config(**kwargs):

    config_def = {
        'creator': creator.base.SaveToCommon,
        'order': ['input', 'common', 'spec_win', 'spectrum_sampling', 'instrument', 'atmosphere', 'radiative_transfer', 'forward_model' , 'retrieval'],
        'input': {
            'creator': creator.base.SaveToCommon,
            'l1b': rf.ExampleLevel1b(l1b_file, observation_id),
            'met': rf.ExampleMetFile(met_file, observation_id),
        },
        'common': {
            'creator': creator.base.SaveToCommon,
            'desc_band_name': static_value("Common/desc_band_name", dtype=str),
            'hdf_band_name': static_value("Common/hdf_band_name", dtype=str),
            'band_reference': {
                'creator': creator.value.ArrayWithUnit,
                'value': static_value("Common/band_reference_point"),
                'units': static_units("Common/band_reference_point"),
            },
            'num_channels': num_channels,
            'absco_base_path': os.environ['ABSCO_PATH'],
            'constants': {
                'creator': creator.common.DefaultConstants,
            },
            'stokes_coefficient': {
                'creator': creator.l1b.ValueFromLevel1b,
                'field': "stokes_coefficient",
            },
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
            'creator': creator.forward_model.FixedSpacingSpectrumSampling,
            'high_res_spacing': rf.DoubleWithUnit(0.01, "cm^-1"), 
        },
        'instrument': {
            'creator': creator.instrument.IlsGratingInstrument,
            'ils_half_width': {
                'creator': creator.value.ArrayWithUnit,
                'value': np.array([4.09e-04, 1.08e-03, 1.40e-03]),
                'units': "um",
            },
            'dispersion': {
                'creator': creator.instrument.DispersionPolynomial,
                'polynomial_coeffs': {
                    'creator': creator.l1b.ValueFromLevel1b,
                    'field': 'spectral_coefficient',
                },
                'number_samples': static_value("Instrument/Dispersion/number_pixel"),
                'is_one_based': True,
                'num_parameters': 2,
            },
            'ils_function': {
                'creator': creator.instrument.IlsTable,
                'delta_lambda': ils_delta_lambda(),
                'response': ils_relative_response(),
            },
            'instrument_correction': {
                'creator': creator.instrument.InstrumentCorrectionList,
                'corrections': [],
            },
        },
        'atmosphere': {
            'creator': creator.atmosphere.AtmosphereCreator,
            'pressure': {
                'creator': creator.atmosphere.PressureSigma,
                'surface_pressure': {
                    'creator': creator.met.ValueFromMet,
                    'field': "surface_pressure",
                },
                'a_coeff': static_value("Pressure/Pressure_sigma_a"),
                'b_coeff': static_value("Pressure/Pressure_sigma_b"),
            },
            'temperature': {
                'creator': creator.atmosphere.TemperatureMet,
                'offset': static_scalar("Temperature/Offset/a_priori")
            },
            'altitude': { 
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
                'creator': creator.absorber.AbsorberAbsco,
                'gases': ['CO2', 'H2O', 'O2'],
                'CO2': {
                    'creator': creator.absorber.AbsorberGasDefinition,
                    'vmr': {
                        'creator': creator.absorber.AbsorberVmrLevel,
                        'vmr_profile': {
                            'creator': creator.absorber.GasVmrAprioriMetL1b,
                            'reference_atm_file': static_input_file,
                        },
                    },
                    'absorption': {
                        'creator': creator.absorber.AbscoLegacy,
                        'table_scale': [1.0, 1.0, 1.004],
                        'filename': "{absco_base_path}/co2_devi2015_wco2scale-nist_sco2scale-unity.h5",
                    },
                },
                'H2O': {
                    'creator': creator.absorber.AbsorberGasDefinition,
                    'vmr': {
                        'creator': creator.absorber.AbsorberVmrMet,
                        'vmr_profile': np.array([1.0]),
                    },
                    'absorption': {
                        'creator': creator.absorber.AbscoLegacy,
                        'table_scale': 1.0,
                        'filename': "{absco_base_path}/h2o_hitran12.h5",
                    },
                },
                'O2': {
                    'creator': creator.absorber.AbsorberGasDefinition,
                    'vmr': {
                        'creator': creator.absorber.AbsorberVmrLevel,
                        'vmr_profile': {
                            'creator': creator.atmosphere.ConstantForAllLevels,
                            'value': static_value("Gas/O2/average_mole_fraction")[0],
                        },
                        'retrieved': False,
                    },
                    'absorption': {
                        'creator': creator.absorber.AbscoLegacy,
                        'table_scale': 1.0,
                        'filename': "{absco_base_path}/o2_v151005_cia_mlawer_v151005r1_narrow.h5",
                     },
                },
            },
            'aerosol': {
                'creator': creator.aerosol.AerosolOptical,
                'aerosols': [ "kahn_2b", "kahn_3b", "water", "ice" ],
                'kahn_2b': {
                    'creator': creator.aerosol.AerosolDefinition,
                    'extinction': {
                        'creator': creator.aerosol.AerosolShapeGaussian,
                        'shape_params': np.array([-4.38203, 1, 0.2]),
                    },
                    'properties': {
                        'creator': creator.aerosol.AerosolPropertyHdf,
                        'filename': aerosol_prop_file,
                    },
                },
                'kahn_3b': {
                    'creator': creator.aerosol.AerosolDefinition,
                    'extinction': {
                        'creator': creator.aerosol.AerosolShapeGaussian,
                        'shape_params': np.array([-4.38203, 1, 0.2]),
                    },
                    'properties': {
                        'creator': creator.aerosol.AerosolPropertyHdf,
                        'filename': aerosol_prop_file,
                    },
                },
                'water': {
                    'creator': creator.aerosol.AerosolDefinition,
                    'extinction': {
                        'creator': creator.aerosol.AerosolShapeGaussian,
                        'shape_params': np.array([-4.38203, 0.75, 0.1]),
                    },
                    'properties': {
                        'creator': creator.aerosol.AerosolPropertyHdf,
                        'filename': aerosol_prop_file,
                        'prop_name': "wc_008",
                    },
                },
                'ice': {
                    'creator': creator.aerosol.AerosolDefinition,
                    'extinction': {
                        'creator': creator.aerosol.AerosolShapeGaussian,
                        'shape_params': np.array([-4.38203, 0.3, 0.04]),
                    },
                    'properties': {
                        'creator': creator.aerosol.AerosolPropertyHdf,
                        'filename': aerosol_prop_file,
                        'prop_name': "ice_cloud_MODIS6_deltaM_1000",
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
                    'polynomial_coeffs': {
                        'creator': creator.ground.AlbedoFromSignalLevel,
                        'signal_level': {
                            'creator': creator.l1b.SignalLevelFromL1b,
                            'sample_grid_indexes': SampleGridIndexes,
                        },
                        'solar_zenith': {
                            'creator': creator.l1b.ValueFromLevel1b,
                            'field': "solar_zenith",
                        },
                        'solar_strength': np.array([4.87e21, 2.096e21, 1.15e21]),
                        'solar_distance': {
                            'creator': creator.l1b.SolarDistanceFromL1b,
                        },
                    },
                },
            },
        },
        'radiative_transfer': {
            'creator': creator.rt.LsiRt,
            'solar_zenith': {
                'creator': creator.l1b.ValueFromLevel1b,
                'field': "solar_zenith",
            },
            'observation_zenith': {
                'creator': creator.l1b.ValueFromLevel1b,
                'field': "sounding_zenith",
            },
            'relative_azimuth': {
                'creator': creator.l1b.ValueFromLevel1b,
                'field': "relative_azimuth",
            },
            'num_low_streams': 1,
            'num_high_streams': 8,
            'lsi_config_file': static_input_file,
        },
        'forward_model': {
            'creator': creator.forward_model.ForwardModel,
            'spectrum_effect': {
                'creator': creator.forward_model.SpectrumEffectList,
                'effects': ["solar_model","instrument_doppler"],
                'solar_model': {
                    'creator': creator.solar_model.SolarAbsorptionAndContinuum,
                    'doppler': {
                        'creator': creator.solar_model.SolarDopplerShiftPolynomial,
                        'time': {
                            'creator': creator.l1b.ValueFromLevel1b,
                            'field': "time",
                        },
                        'latitude': {
                            'creator': creator.l1b.ValueFromLevel1b,
                            'field': "latitude",
                        },
                        'solar_zenith': {
                            'creator': creator.l1b.ValueFromLevel1b,
                            'field': "solar_zenith",
                        },
                        'solar_azimuth': {
                            'creator': creator.l1b.ValueFromLevel1b,
                            'field': "solar_azimuth",
                        },
                        'surface_height': {
                            'creator': creator.l1b.ValueFromLevel1b,
                            'field': "altitude",
                        },
                    },
                    'absorption': {
                        'creator': creator.solar_model.SolarAbsorptionTable,
                        'solar_data_file': solar_file,
                    },
                    'continuum': {
                        'creator': creator.solar_model.SolarContinuumTable,
                        'solar_data_file': solar_file,
                    },
                },
                'instrument_doppler': {
                    'creator': creator.instrument.InstrumentDoppler,
                    'relative_velocity': {
                        'creator': creator.l1b.ValueFromLevel1b,
                        'field': "relative_velocity",
                    },
                },
            },
        },
        'retrieval': {
            'creator': creator.retrieval.NLLSRetrieval,
            'retrieval_components': {
                'creator': creator.retrieval.SVObserverComponents,
                'exclude': ['absorber_levels/linear/O2', 'instrument_doppler'],
                # Match order tradtionally used in old system
                'order': ['CO2', 'H2O', 'surface_pressure', 'temperature_offset', 'aerosol', 'ground', 'dispersion'],
            },
            'state_vector': {
                'creator': creator.retrieval.StateVector,
            },
            'initial_guess': {
                'creator': creator.retrieval.InitialGuessFromSV,
            },
            'a_priori': {
                'creator': creator.retrieval.AprioriFromIG,
            },
            'covariance': {
                'creator': creator.retrieval.CovarianceByComponent,
                'values': {
                    'creator': creator.value.LoadValuesFromHDF,
                    'filename': covariance_file,
                }
            },
            'solver': {
                'creator': creator.retrieval.LegacyConnorSolver,
                'max_iteration': 7,
                'max_divergence': 2,
                'max_chisq': 1.4,
                'threshold': 2.0,
                'gamma_initial': 10.0,
            },
        },
    }

    return config_def
