#ifndef OSS_CONFIGURATION_FIXTURE_H
#define OSS_CONFIGURATION_FIXTURE_H
/*
#include "absorber.h"
#include "aerosol_optical.h"

#include "error_analysis.h"
#include "connor_solver.h"
#include "spectral_window.h"
#include "spectrum_sampling.h"
#include "level_1b_sample_coefficient.h"
#include "initial_guess.h"
#include "instrument.h"
#include "sample_grid.h"
*/
#include <blitz/array.h>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

#include "global_fixture.h"
#include "pressure_level_input.h"
#include "pressure_fixed_level.h"
#include "temperature_fixed_level.h"
#include "rt_atmosphere.h"
#include "hdf_file.h"
#include "absorber.h"
/*
#include "radiative_transfer.h"
#include "lua_state.h"

#include "forward_model.h"
#include "observation.h"
#include <map>
*/

namespace FullPhysics {
/****************************************************************//**
  The OSS FM uses a different test atmosphere than the rest of the
  framework. This file is only used by OSS FM tests.
*******************************************************************/
class OssConfigurationFixture: public GlobalFixture {
public:

  /* Note: input_file is expected to be within oss_data_dir() */
  OssConfigurationFixture(const std::string& input_file = "tape5_nc4.nc");
  ~OssConfigurationFixture();
  /*
  boost::shared_ptr<Absorber> config_absorber;
  boost::shared_ptr<Aerosol> config_aerosol;
  */
  boost::shared_ptr<Absorber> config_absorber;
  std::vector<bool> config_absorber_calc_jacob;
  boost::shared_ptr<Pressure> config_pressure;
  boost::shared_ptr<Temperature> config_temperature;
  DoubleWithUnit config_skin_temperature;
  boost::shared_ptr<RtAtmosphere> config_atmosphere;

  /*
  boost::shared_ptr<StateVector> config_state_vector;
  boost::shared_ptr<PressureLevelInput> config_pressure_level_input;
  boost::shared_ptr<ErrorAnalysis> config_error_analysis;
  boost::shared_ptr<ConnorSolver> config_solver;
  boost::shared_ptr<SpectralWindow> config_spectral_window;
  boost::shared_ptr<InitialGuess> config_initial_guess;
  boost::shared_ptr<Instrument> config_instrument;


  boost::shared_ptr<SpectrumSampling> config_spectrum_sampling;
  boost::shared_ptr<Level1bSampleCoefficient> config_level_1b;
  boost::shared_ptr<Ground> config_ground;
  boost::shared_ptr<ForwardModel> config_forward_model;
  boost::shared_ptr<Observation> config_observation;
  boost::shared_ptr<RadiativeTransfer> config_rt;
  */


private:
  boost::shared_ptr<HdfFile> input_data;
};
}
#endif
