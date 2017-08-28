
-- This sets up the "standard" run we use in the unit tests, you can then
-- use this and override any feature you want to be done differently

require "config_common"
base_config_dir = ConfigCommon.local_dir()

local example_level1b = Creator:new()

function example_level1b:create()
   local hv = self.config:l1b_hdf_file()
   return ExampleLevel1b(hv, self.config.sid_string)
end

local example_met = Creator:new()

function example_met:create()
    local met_hdf = HdfFile(self.config.met_file)
    return ExampleMetFile(met_hdf, self.config.sid_string)
end

local ils_table_creator = Creator:new()

function ils_table_creator:create()
    -- Use a simple table ILS since the ILS supplies
    -- a full set of delta_lambda/repsonses per pixel
    interpolate = false

    local ils_hdf_file = HdfFile(self.config.static_ils_data_file)

    local idx
    local wavenumber, delta_lambda, response_name
    local res = {}
    for idx = 0, self.config.number_pixel:rows() - 1 do
        local hdf_band_name = self.config.common.hdf_band_name:value(idx)
        local desc_band_name = self.config.common.desc_band_name:value(idx)

        -- Delta lambda
        delta_lambda = ils_hdf_file:read_double_3d("/InstrumentData/ils_delta_lambda")
        delta_lambda = delta_lambda(idx, Range.all(), Range.all())

        -- Load ILS response values
        response = ils_hdf_file:read_double_3d("/InstrumentData/ils_relative_response")
        response = response(idx, Range.all(), Range.all())

        -- Calculate center wavenumber from dispersion, there should be number pixel
        -- of these per spectrometer
        wavenumber = self.config.dispersion[idx+1]:pixel_grid():data()

        res[idx+1] = IlsTableLog(wavenumber, delta_lambda, response, desc_band_name, hdf_band_name, interpolate)
    end
    return res
end

ExampleBaseConfig = ConfigCommon:new {
------------------------------------------------------------
--- Various constants used to describe input data.
------------------------------------------------------------

    sid_string = "2014090915251774",
    spectrum_file = base_config_dir .. "/../in/common/l1b_example_data.h5",
    met_file = base_config_dir .. "/../in/common/met_example_data.h5",

    static_file = base_config_dir .. "/example_static_input.h5",
    static_ils_data_file = base_config_dir .. "/ils_data.h5",

    static_solar_file = config_common_dir .. "/../input/l2_solar_model.h5",
    static_aerosol_file = config_common_dir .. "/../input/l2_aerosol_combined.h5",

------------------------------------------------------------
-- Set this true to get diagnostic messages to help debug
-- problems with Lua
------------------------------------------------------------

    diagnostic = false,

------------------------------------------------------------
-- Paths for Absco data. We first look in local path, and if
-- we can't read the file try the full path. This allows us
-- to use a local disk on the cluster.
------------------------------------------------------------

    absco_path = "/path/to/absco",

------------------------------------------------------------
-- Connor solver
------------------------------------------------------------

    solver = { 
        max_iteration=7, max_divergence=2,
        max_chisq=1.4, threshold=2.0, gamma_initial=10.0,
        create = ConfigCommon.connor_solver,
    },

------------------------------------------------------------
-- If true then launch solver, otherwise just do a 
-- forward model calculation, with jacobians if
-- write_jacobians is true
------------------------------------------------------------

    do_retrieval = true,

------------------------------------------------------------
-- True if we want to write the jacobians out
------------------------------------------------------------

    write_jacobian = false,

------------------------------------------------------------
-- True if we want to write high resolution spectra out
------------------------------------------------------------

    write_high_res_spectra = false,

------------------------------------------------------------
-- True if we want to generate output for every iteration
------------------------------------------------------------

    iteration_output = false,

------------------------------------------------------------
--- Log level
------------------------------------------------------------

    log_level = LogImp.INFO,

-----------------------------------------------------------
--- Default creators for everything. We default to getting stuff
--- from the HDF file, but this can be changed if desired.
------------------------------------------------------------

    fm = {
        creator = ConfigCommon.oco_forward_model,
        common = {
            desc_band_name = ConfigCommon.hdf_read_string_vector("Common/desc_band_name"),
            hdf_band_name = ConfigCommon.hdf_read_string_vector("Common/hdf_band_name"),
            band_reference = ConfigCommon.hdf_read_double_with_unit_1d("Common/band_reference_point"),
            creator = ConfigCommon.table_function_eval,
        },
        spec_win = {
            creator = ConfigCommon.spectral_window_hdf,
        },
        input = {
            creator = ConfigCommon.l1b_met_input,
            l1b = {
                creator = example_level1b,
            },
            met = {
                creator = example_met,
            },
        },
        stokes_coefficient = {
            creator = ConfigCommon.stokes_coefficient_constant,
            value = ConfigCommon.stokes_coefficient_l1b,
        },
        instrument = {
            creator = ConfigCommon.ils_instrument,
            ils_half_width = { DoubleWithUnit(4.09e-04, "um"), 
                               DoubleWithUnit(1.08e-03, "um"),
                               DoubleWithUnit(1.40e-03, "um") },
            dispersion = {
                creator = ConfigCommon.dispersion_polynomial,
                apriori = ConfigCommon.l1b_spectral_coefficient_i,
                covariance = ConfigCommon.hdf_covariance_i("Instrument/Dispersion"),
                number_pixel = ConfigCommon.hdf_read_int_1d("Instrument/Dispersion/number_pixel"),
                retrieved = true,
                is_one_based = true,
                num_parameters = 2,
            },
            ils_func = {
                creator = ils_table_creator,
            },
            instrument_correction = {
                creator = ConfigCommon.instrument_correction_list,
                ic = {},
            },
        },
        spectrum_effect = {
            creator = ConfigCommon.spectrum_effect_list,
            speceff = { "solar_model", "instrument_doppler" },
            solar_model = {
                creator = ConfigCommon.solar_absorption_and_continuum,
                doppler_shift = {
                    creator = ConfigCommon.solar_doppler_from_l1b,
                    do_doppler_shift = true,
                },
                solar_absorption = {
                    creator = ConfigCommon.solar_absorption_table,
                },
                solar_continuum = {
                    creator = ConfigCommon.solar_continuum_table,
                    convert_from_photon = false,
                },
            },
            instrument_doppler = {
                creator = ConfigCommon.instrument_doppler,
                retrieved = false,
            },
        },
        spec_samp = {
            creator = ConfigCommon.uniform_spectrum_sampling,
            high_resolution_spectrum_spacing = DoubleWithUnit(0.01, "cm^-1"),
        },
        rt = {
            creator = ConfigCommon.radiative_transfer_lsi,
            nadir_threshold = 1e-6,
            lsi_constant = {
                dedicated_twostream = true,
                -- Note that the "1" here is just a convention to use the
                -- dedicated two stream code
                low_stream = 1,

                -- LIDORT input is in Half-Streams. Full-streams is double
                -- this (so high_stream = 8 would mean 16 full-streams)
                high_stream = 8
            },
        },
        state_vector = {
            creator = ConfigCommon.state_vector_creator,
        },
        atmosphere = {
            creator = ConfigCommon.atmosphere_oco,
            constants = {
                creator = ConfigCommon.default_constant,
            },
            pressure = {
                apriori = ConfigCommon.met_pressure,
                covariance = ConfigCommon.hdf_covariance("Surface_Pressure"),
                a = ConfigCommon.hdf_read_double_1d("Pressure/Pressure_sigma_a"),
                b = ConfigCommon.hdf_read_double_1d("Pressure/Pressure_sigma_b"),
                creator = ConfigCommon.pressure_sigma,
            },
            temperature = {
                apriori = ConfigCommon.hdf_apriori("Temperature/Offset"),
                covariance = ConfigCommon.hdf_covariance("Temperature/Offset"),
                creator = ConfigCommon.temperature_met,
            },
            ground = {
                -- Instrument specific solar strengths used for ground calculations 
                solar_strength = {4.87e21, 2.096e21, 1.15e21},

                -- Pure lambertian
                lambertian = {
                    apriori = ConfigCommon.albedo_from_signal_level(1),
                    covariance = ConfigCommon.hdf_covariance_i("Ground/Albedo"),
                    retrieve_bands = { true, true, true },
                    creator = ConfigCommon.lambertian_retrieval,
                },

                -- Coxmunk windspeed and refractive index inputs
                coxmunk = {
                    refractive_index = ConfigCommon.hdf_apriori("Ground/Refractive_Index"),
                    apriori = ConfigCommon.met_windspeed,
                    covariance = ConfigCommon.hdf_covariance("Ground/Windspeed"),
                    creator = ConfigCommon.coxmunk_retrieval,
                    },

                -- Lambertian component of coxmunk + lambertian
                coxmunk_lambertian = {
                    apriori = ConfigCommon.hdf_apriori_i("Ground/Coxmunk_Albedo"),
                    covariance = ConfigCommon.hdf_covariance_i("Ground/Coxmunk_Albedo"),
                    retrieve_bands = { true, true, true },
                    creator = ConfigCommon.lambertian_retrieval,
                },

                -- Brdf vegetative kernel with Rahman retrieved parameters
                brdf_veg = {
                    apriori = ConfigCommon.brdf_veg_apriori("Ground/Brdf"),
                    covariance = ConfigCommon.hdf_covariance_i("Ground/Brdf"),
                    retrieve_bands = { true, true, true },
                    creator = ConfigCommon.brdf_veg_retrieval,
                },
            
                -- Brdf soil kernel with Rahman retrieved parameters
                brdf_soil = {
                    apriori = ConfigCommon.brdf_soil_apriori("Ground/Brdf"),
                    covariance = ConfigCommon.hdf_covariance_i("Ground/Brdf"),
                    retrieve_bands = { true, true, true },
                    creator = ConfigCommon.brdf_soil_retrieval,
                },

                creator = ConfigCommon.ground_lambertian,
            },
            aerosol = {
                creator = ConfigCommon.aerosol_creator,
                -- Lua doesn't preserve order in a table, so we have a list
                -- saying what order we want the Aerosols in
                aerosols = {"Kahn_2b", "Kahn_3b", "Water", "Ice"},
                Kahn_2b = {
                    creator = ConfigCommon.aerosol_log_shape_gaussian,
                    apriori = ConfigCommon.hdf_aerosol_apriori("Aerosol", "Gaussian/Log"),
                    covariance = ConfigCommon.hdf_aerosol_covariance("Aerosol", "Gaussian/Log"),
                    property = ConfigCommon.hdf_aerosol_property("kahn_2b"),
                },
                Kahn_3b = {
                    creator = ConfigCommon.aerosol_log_shape_gaussian,
                    apriori = ConfigCommon.hdf_aerosol_apriori("Aerosol", "Gaussian/Log"),
                    covariance = ConfigCommon.hdf_aerosol_covariance("Aerosol", "Gaussian/Log"),
                    property = ConfigCommon.hdf_aerosol_property("kahn_3b"),
                },
                Water = {
                    creator = ConfigCommon.aerosol_log_shape_gaussian,
                    apriori = ConfigCommon.hdf_aerosol_apriori("Aerosol", "Gaussian/Log"),
                    covariance = ConfigCommon.hdf_aerosol_covariance("Aerosol", "Gaussian/Log"),
                    property = ConfigCommon.hdf_aerosol_property("wc_008"),
                },
                Ice = {
                    creator = ConfigCommon.aerosol_log_shape_gaussian,
                    apriori = ConfigCommon.hdf_aerosol_apriori("Aerosol", "Gaussian/Log"),
                    covariance = ConfigCommon.hdf_aerosol_covariance("Aerosol", "Gaussian/Log"),
                    property = ConfigCommon.hdf_aerosol_property("ice_cloud_MODIS6_deltaM_1000"),
                },
            },
            absorber = {
                creator = ConfigCommon.absorber_creator,
                gases = {"CO2", "H2O", "O2"},
                CO2 = {
                    apriori = ConfigCommon.reference_co2_apriori_met_apriori,
                    covariance = ConfigCommon.hdf_covariance("Gas/CO2"),
                    absco = "v5.0.0/co2_devi2015_wco2scale-nist_sco2scale-unity.h5",
                    table_scale = {1.0, 1.0, 1.004},
                    creator = ConfigCommon.vmr_level,
                },
                H2O = {
                    scale_apriori = 1.0,
                    scale_cov = 0.25,
                    absco = "v5.0.0/h2o_hitran12.h5",
                    creator = ConfigCommon.vmr_met,
                },
                O2 = {
                    apriori = ConfigCommon.hdf_read_double_1d("Gas/O2/average_mole_fraction"),
                    absco = "v5.0.0/o2_v151005_cia_mlawer_v151005r1_narrow.h5",
                    table_scale = 1.0,
                    creator = ConfigCommon.vmr_level_constant_well_mixed,
                },
            },
            altitude = {
                creator = ConfigCommon.hydrostatic_altitude,
            },
            relative_humidity = {
                creator = ConfigCommon.calc_relative_humidity,
            },
        }, -- end of atmosphere
    }, -- end of fm
} -- end of ExampleBaseConfig
