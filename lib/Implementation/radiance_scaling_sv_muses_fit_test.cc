#include "unit_test_support.h"

#include "serialized_configuration_fixture.h"
#include "radiance_scaling_sv_muses_fit.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(radiance_scaling_sv_muses_fit, LambertianConfigurationFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  blitz::Array<double, 1> coeff(2);
  coeff = 0.5, 0.0;
  DoubleWithUnit band_ref(765, "nm");
  std::string band_name = "ABand";
  RadianceScalingSvMusesFit rs1(coeff, band_ref, band_name);
  
  SpectralDomain pixel_grid = config_instrument->pixel_spectral_domain(0);
  std::vector<int> pixel_list;
  for(int i = 0; i < pixel_grid.data().rows(); i++)
    pixel_list.push_back(i);
  
  blitz::Array<double, 1> radiance_data(pixel_grid.data().shape());
  radiance_data = 1.0;
  SpectralRange radiance(radiance_data, Unit());

  rs1.apply_correction(pixel_grid, pixel_list, radiance);

  blitz::Array<double, 1> rad_expect(radiance_data.shape());
  rad_expect = 0.5;
  BOOST_CHECK_MATRIX_CLOSE(radiance.data(), rad_expect);

  // Check array ad and more complex coeff
  ArrayAd<double, 1> radiance_data2(pixel_grid.data().shape(), 5);
  radiance_data2.value() = 1.0;
  radiance_data2.jacobian() = 0.0;
  SpectralRange radiance2(radiance_data2, Unit());

  coeff = 0.5, 1e-3;
  RadianceScalingSvMusesFit rs2(coeff, band_ref, band_name);

  // Set up statevector stuff so that we can properly
  // test the jacobians going through radiance scaling
  StateVector sv;
  sv.add_observer(rs2);
  Array<double,1> x(5);
  x(Range(0,1)) = coeff(Range::all()); 
  sv.update_state(x);

  rs2.apply_correction(pixel_grid, pixel_list, radiance2);

  double band_ref_conv = band_ref.convert_wave(Unit("nm")).value;
  Array<double, 1> gd = pixel_grid.convert_wave(Unit("nm"));
  
  for(int i = 0; i < pixel_grid.data().rows(); i++)
    rad_expect(i) = coeff(0) + coeff(1) * (1 - gd(i) / band_ref_conv);
  for(int i = 0; i < pixel_grid.data().rows(); i++)
  BOOST_CHECK_MATRIX_CLOSE_TOL(radiance2.data(), rad_expect, 1e-8);
  
  for(int i = 0; i < pixel_grid.data().rows(); i++) {
    BOOST_CHECK_CLOSE(radiance2.data_ad().jacobian()(i, 0), 1, 1e-7);
    BOOST_CHECK_CLOSE(radiance2.data_ad().jacobian()(i, 1), 1 - gd(i)/ band_ref_conv, 1e-7);
  }
}

BOOST_AUTO_TEST_CASE(serialization)
{
  if(!have_serialize_supported())
    return;
  blitz::Array<double, 1> coeff(2);
  coeff = 0.5, 0.0;
  DoubleWithUnit band_ref(0.77, "cm^-1");
  std::string band_name = "ABand";
  boost::shared_ptr<RadianceScalingSvMusesFit> rs = boost::make_shared<RadianceScalingSvMusesFit>(coeff, band_ref, band_name);
  std::string d = serialize_write_string(rs);
  if(false)
    std::cerr << d;
  boost::shared_ptr<RadianceScalingSvMusesFit> rsr = serialize_read_string<RadianceScalingSvMusesFit>(d);
  SpectralDomain pixel_grid = config_instrument->pixel_spectral_domain(0);
  std::vector<int> pixel_list;
  for(int i = 0; i < pixel_grid.data().rows(); i++)
    pixel_list.push_back(i);
  
  blitz::Array<double, 1> radiance_data(pixel_grid.data().shape());
  radiance_data = 1.0;
  SpectralRange radiance(radiance_data, Unit());

  rsr->apply_correction(pixel_grid, pixel_list, radiance);

  blitz::Array<double, 1> rad_expect(radiance_data.shape());
  rad_expect = 0.5;
  BOOST_CHECK_MATRIX_CLOSE(radiance.data(), rad_expect);
}

BOOST_AUTO_TEST_SUITE_END()
