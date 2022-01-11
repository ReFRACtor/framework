#include "serialized_configuration_fixture.h"

// Unsure of why but without including unit_test_support.h the import of p_serialize_support.h fails with
// complaints about GenericObject not being defined?!
#include "unit_test_support.h"
#include "generic_object_map.h"
#include "fp_serialize_support.h"

using namespace blitz;
using namespace FullPhysics;

SerializedConfigurationFixture::SerializedConfigurationFixture(const std::string& Serialized_file)
  
{
    // Preprend test_data_dir here this function is only available from within the scope of a fixture
    serialized_filename = test_data_dir() + "in/configuration_fixture/" + Serialized_file;

    init_variables();
    init_epsilon();
}

void SerializedConfigurationFixture::init_variables()
{
    boost::shared_ptr<GenericObjectMap> obj_map = serialize_read_binary<GenericObjectMap>(serialized_filename);

    config_absorber = obj_map->get<Absorber>("absorber");
    config_pressure = obj_map->get<Pressure>("pressure");
    config_aerosol = obj_map->get<Aerosol>("aerosol");
    config_atmosphere = obj_map->get<RtAtmosphere>("atmosphere");
    config_state_vector = obj_map->get<StateVector>("state_vector");
    config_solver = obj_map->get<ConnorSolver>("solver");
    config_spectral_window = obj_map->get<SpectralWindow>("spectral_window");
    config_initial_guess = obj_map->get<InitialGuess>("initial_guess");
    config_instrument = obj_map->get<Instrument>("instrument");
    config_temperature = obj_map->get<Temperature>("temperature");
    config_spectrum_sampling = obj_map->get<SpectrumSampling>("spectrum_sampling");
    config_level_1b = obj_map->get<Level1bSampleCoefficient>("level_1b");
    config_ground = obj_map->get<Ground>("ground");
    config_forward_model = obj_map->get<ForwardModel>("forward_model");
    config_observation = obj_map->get<Observation>("observation");
    config_rt = obj_map->get<RadiativeTransfer>("rt");

}

void SerializedConfigurationFixture::init_epsilon()
{
    epsilon.resize(config_state_vector->observer_claimed_size());
    epsilon = 1e-6;                  // Default

    for (int sv_idx = 0 ; sv_idx < config_state_vector->observer_claimed_size(); sv_idx++) {
        std::string sv_name = config_state_vector->state_vector_name()(sv_idx);

        if (sv_name.find("CO2") == 0) {
            // CO2 VMR
            epsilon(sv_idx) = 1e-7;
        } else if (sv_name.find("Surface Pressure") == 0) {
            // Surface Pressure
            epsilon(sv_idx) = 1e-3;
        } else if (sv_name.find("Temperature") == 0) {
            // Temperature
            epsilon(sv_idx) = 1e-4;
        } else if (sv_name.find("Aerosol") == 0) {
            // Aerosol
            epsilon(sv_idx) = 1e-8;
        }

    }
}
