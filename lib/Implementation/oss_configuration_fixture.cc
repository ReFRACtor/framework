#include "oss_configuration_fixture.h"

using namespace FullPhysics;
using namespace blitz;


OssConfigurationFixture::OssConfigurationFixture(const std::string& input_file)
{
  input_data = boost::make_shared<HdfFile>(oss_run_dir() + input_file);

  /* Pressure in file from high to low */
  Array<float, 1> pressure_nc = input_data->read_field<float, 1>("/Pressure")(Range::all());

  double surface_pressure = static_cast<double>(pressure_nc(0));

  Array<double, 1> pressure_above_surface(pressure_nc.rows() - 1);
  for (int i = 1; i < pressure_nc.rows(); i++) {
    pressure_above_surface(i - 1) = static_cast<double>(pressure_nc(i));
  }

  boost::shared_ptr<PressureLevelInput> pressure_level = boost::make_shared<PressureLevelInput>(pressure_above_surface);

  config_pressure = boost::make_shared<PressureFixedLevel>(true, pressure_level, surface_pressure);


  /* Temperature */
  Array<float, 1> temp_nc = input_data->read_field<float, 1>("/Temperature")(Range::all());

  double surface_temp = static_cast<double>(temp_nc(0));

  Array<double, 1> temp(temp_nc.rows() - 1);
  for (int i = 1; i < temp_nc.rows(); i++) {
    temp(i - 1) = static_cast<double>(temp_nc(i));
  }
  Array<bool, 1> temp_flag(temp_nc.rows() - 1);
  temp_flag = true;
/*
 TemperatureFixedLevel::TemperatureFixedLevel
(const blitz::Array<bool, 1>& Flag_temp,
 bool Flag_offset, const blitz::Array<double, 1>& Temp,
 double T_offset,
 const boost::shared_ptr<Pressure>& Press,
 const boost::shared_ptr<PressureLevelInput>& Press_level)

 */
  config_temperature = boost::make_shared<TemperatureFixedLevel>(temp_flag, true, temp, surface_temp, config_pressure, pressure_level);


  /*
  AtmosphereStandard(const boost::shared_ptr<Absorber>& absorberv,
                     const boost::shared_ptr<Pressure>& pressurev,
                     const boost::shared_ptr<Temperature>& temperaturev,
                     const boost::shared_ptr<Aerosol>& aerosolv,
                     const boost::shared_ptr<RelativeHumidity>& rhv,
                     const boost::shared_ptr<Ground>& groundv,
                     const boost::shared_ptr<SurfaceTemperature>& surface_tempv,
                     const std::vector<boost::shared_ptr<Altitude> >& altv,
                     const boost::shared_ptr<Constant>& C);
   */
}

OssConfigurationFixture::~OssConfigurationFixture() {};
