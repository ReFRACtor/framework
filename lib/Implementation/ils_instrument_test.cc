#include "ils_instrument.h"
#include "ils_grating.h"
#include "hdf_file.h"
#include "dispersion_polynomial.h"
#include "sample_grid_spectral_domain.h"
#include "ils_table.h"
#include "unit_test_support.h"
#include "serialized_configuration_fixture.h"
#include "accumulated_timer.h"
#include <iostream>
using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(ils_instrument, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  HdfFile hf(test_data_dir() + "in/ils/ils_linear_table.h5");
  Array<double, 1> disp_coeff(2);
  boost::shared_ptr<IlsTableLinear> ils_tab(new IlsTableLinear(hf, 0, "A-Band", "o2"));

  disp_coeff = 1.28695614e+04, 1.99492886e-01;
  firstIndex i1;
  Array<double, 1> disp_var_1(1805);
  disp_var_1 = i1 + 1;

  boost::shared_ptr<DispersionPolynomial>
    d(new DispersionPolynomial(disp_coeff, units::inv_cm, disp_var_1, ils_tab->band_name()));

  std::vector<boost::shared_ptr<Ils> > ils;
  ils.push_back(boost::shared_ptr<Ils>(new IlsGrating(d, ils_tab)));
  
  ils_tab.reset(new IlsTableLinear(hf, 1, "WC-Band", "weak_co2"));

  disp_coeff = 5.74982835e+03, 1.99492886e-01;
  Array<double, 1> disp_var_2(3508);
  disp_var_2 = i1 + 1;

  d.reset(new DispersionPolynomial(disp_coeff, units::inv_cm, disp_var_2, ils_tab->band_name()));
  ils.push_back(boost::shared_ptr<Ils>(new IlsGrating(d, ils_tab)));

  ils_tab.reset(new IlsTableLinear(hf, 2, "SC-Band", "strong_co2"));

  disp_coeff = 4.74980283e+03, 1.99492886e-01;
  Array<double, 1> disp_var_3(2005);
  disp_var_3 = i1 + 1;

  d.reset(new DispersionPolynomial(disp_coeff, units::inv_cm, disp_var_3, ils_tab->band_name()));

  ils.push_back(boost::shared_ptr<Ils>(new IlsGrating(d, ils_tab)));

  IlsInstrument inst(ils);
  StateVector sv;
  sv.add_observer(inst);
  Array<double,1> x(4);
  x = 1.28695614e+04, 5.74982835e+03,4.74980283e+03,0;
  sv.update_state(x);

  std::vector<int> plist;
  // Arbitrary list of pixels.
  plist.push_back(403);
  plist.push_back(405);
  for(int i = 407; i <= 414; ++i)
    plist.push_back(i);

  IfstreamCs expected(test_data_dir() + "expected/ils_grating/basic");
  Array<double, 1> wn_in, rad_hres_in, rad_out_expect;
  expected >> wn_in >> rad_hres_in >> rad_out_expect;
  SpectralDomain sd(wn_in, units::inv_cm);
  SpectralRange sr(rad_hres_in, units::dimensionless); // Ignore units
  BOOST_CHECK_MATRIX_CLOSE(inst.apply_instrument_model(Spectrum(sd, sr),
                                                       plist, 0).
                           spectral_range().data(), rad_out_expect);

  Array<double, 2> jac_rad_fake(rad_hres_in.rows(), 4);
  jac_rad_fake = 0;
  jac_rad_fake(Range::all(), 3) = rad_hres_in;
  ArrayAd<double, 1> rad_hres_in2(rad_hres_in, jac_rad_fake);
  SpectralRange sr2(rad_hres_in2, units::dimensionless); // Ignore units
  BOOST_CHECK_MATRIX_CLOSE(inst.apply_instrument_model(Spectrum(sd, sr), 
                                                       plist, 0).
                           spectral_range().data(), rad_out_expect);
  Array<double, 2> jac = 
    inst.apply_instrument_model(Spectrum(sd, sr2), plist, 0).
    spectral_range().data_ad().jacobian();
  Array<double, 1> v0 = 
    inst.apply_instrument_model(Spectrum(sd, sr2), plist, 0).
    spectral_range().data();
  double epsilon = 1e-3;
  x(0) += epsilon;
  sv.update_state(x);
  Array<double, 1> v1 = 
    inst.apply_instrument_model(Spectrum(sd, sr2), plist, 0).
    spectral_range().data();
  Array<double, 2> jacd(jac.shape());
  jacd = 0;
  jacd(Range::all(), 0) = (v1 - v0) / epsilon;
  x(0) -= epsilon;
  sv.update_state(x);
  rad_hres_in2.value() *= 1 + epsilon;
  SpectralRange sr3(rad_hres_in2, units::dimensionless); // Ignore units
  Array<double, 1> v2 = 
    inst.apply_instrument_model(Spectrum(sd, sr2), plist, 0).
    spectral_range().data();
  jacd(Range::all(), 3) = (v2 - v0) / epsilon;
  BOOST_CHECK_MATRIX_CLOSE(jac, jacd);
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(ils_instrument_static, LambertianConfigurationFixture)

BOOST_AUTO_TEST_CASE(basic_static_test)
{
  HdfFile hf(test_data_dir() + "in/ils/ils_linear_table.h5");
  boost::shared_ptr<IlsTableLinear> ils_tab(new IlsTableLinear(hf, 0, "A-Band", "o2"));
  IfstreamCs sd_in_data(test_data_dir() + "expected/dispersion_polynomial/gosat");
  Array<double, 1> pwn_in;
  sd_in_data >> pwn_in;
  SpectralDomain sd_in(pwn_in, units::inv_cm);
  boost::shared_ptr<SampleGridSpectralDomain>
    sg(new SampleGridSpectralDomain(sd_in, ils_tab->band_name()));
  std::vector<boost::shared_ptr<Ils> > ils;
  ils.push_back(boost::shared_ptr<Ils>(new IlsGrating(sg, ils_tab)));

  IlsInstrument inst(ils);

  std::vector<int> plist;
  // Arbitrary list of pixels.
  plist.push_back(403);
  plist.push_back(405);
  for(int i = 407; i <= 414; ++i)
    plist.push_back(i);

  IfstreamCs expected(test_data_dir() + "expected/ils_grating/basic");
  Array<double, 1> wn_in, rad_hres_in, rad_out_expect;
  expected >> wn_in >> rad_hres_in >> rad_out_expect;
  SpectralDomain sd(wn_in, units::inv_cm);
  SpectralRange sr(rad_hres_in, units::dimensionless); // Ignore units
  BOOST_CHECK_MATRIX_CLOSE(inst.apply_instrument_model(Spectrum(sd, sr),
                                                       plist, 0).
                           spectral_range().data(), rad_out_expect);
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(ils_instrument_diff_static, LambertianConfigurationFixture)

BOOST_AUTO_TEST_CASE(basic_static_different_units_test)
{
  HdfFile hf(test_data_dir() + "in/ils/ils_linear_table.h5");
  boost::shared_ptr<IlsTableLinear> ils_tab(new IlsTableLinear(hf, 0, "A-Band", "o2"));
  IfstreamCs sd_in_data(test_data_dir() + "expected/dispersion_polynomial/gosat");
  Array<double, 1> pwn_in;
  sd_in_data >> pwn_in;
  SpectralDomain sd_in(pwn_in, units::inv_cm);
  boost::shared_ptr<SampleGridSpectralDomain>
    sg(new SampleGridSpectralDomain(sd_in, ils_tab->band_name()));
  std::vector<boost::shared_ptr<Ils> > ils;
  ils.push_back(boost::shared_ptr<Ils>(new IlsGrating(sg, ils_tab)));

  IlsInstrument inst(ils);

  std::vector<int> plist;
  // Arbitrary list of pixels.
  plist.push_back(403);
  plist.push_back(405);
  for(int i = 407; i <= 414; ++i)
    plist.push_back(i);

  IfstreamCs expected(test_data_dir() + "expected/ils_grating/basic");
  Array<double, 1> wn_in, rad_hres_in, rad_out_expect;
  expected >> wn_in >> rad_hres_in >> rad_out_expect;
  SpectralDomain sd(wn_in, units::inv_cm);
  Array<double, 1> sd_convert_data(wn_in.rows());
  sd_convert_data = sd.convert_wave(units::micron).reverse(firstDim)(Range::all());
  SpectralDomain sd_micron(sd_convert_data, units::micron);

  Array<double, 1> sr_convert_data(rad_hres_in.rows());
  sr_convert_data = rad_hres_in.reverse(firstDim)(Range::all());
  SpectralRange sr(sr_convert_data, units::dimensionless); // Ignore units
  BOOST_CHECK_MATRIX_CLOSE(inst.apply_instrument_model(Spectrum(sd_micron, sr),
                                                       plist, 0).
                           spectral_range().data(), rad_out_expect);
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(ils_instrument_timing, LambertianConfigurationFixture)

BOOST_AUTO_TEST_CASE(timing)
{
  is_timing_test();
  const Instrument& inst = *config_instrument;
  IfstreamCs expected(test_data_dir() + "expected/ils_grating/basic");
  Array<double, 1> wn_in, rad_hres_val;
  expected >> wn_in >> rad_hres_val;
  Array<double, 2> rad_hres_jac(rad_hres_val.rows(), 
                                 config_state_vector->state().rows());
  for(int i = 0; i < rad_hres_jac.cols(); ++i)
    rad_hres_jac(Range::all(), i) = rad_hres_val;
  ArrayAd<double, 1> rad_hres(rad_hres_val, rad_hres_jac);
  std::vector<int> plist = 
    config_spectral_window->grid_indexes(inst.pixel_spectral_domain(0), 0);
  Spectrum s(SpectralDomain(wn_in, units::inv_cm),
             SpectralRange(rad_hres, units::dimensionless));
  AccumulatedTimer tm("IlsInstrument");
  {
    FunctionTimer ft = tm.function_timer();
    Spectrum res = inst.apply_instrument_model(s, plist, 0);
  }
  std::cout << tm << "\n";
}

BOOST_AUTO_TEST_SUITE_END()

