#include "empirical_orthogonal_function.h"
#include "ils_instrument.h"
#include "ils_grating.h"
#include "hdf_file.h"
#include "dispersion_polynomial.h"
#include "ils_table.h"
#include "unit_test_support.h"
#include "lua_configuration_fixture.h"
#include "accumulated_timer.h"
#include <iostream>
using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(empirical_orthogonal_function, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  std::vector<std::vector<boost::shared_ptr<InstrumentCorrection> > > corr(3);

  HdfFile hf_ils(test_data_dir() + "in/ils/ils_linear_table.h5");
  HdfFile hf_eof(test_data_dir() + "in/eof/eof.h5");
  Array<double, 1> disp_coeff(2);

  boost::shared_ptr<IlsTableLinear> ils_tab(new IlsTableLinear(hf_ils, 0, "A-Band", "o2"));

  disp_coeff = 1.28695614e+04, 1.99492886e-01;
  firstIndex i1;
  Array<double, 1> disp_var_1(1805);
  disp_var_1 = i1 + 1;

  boost::shared_ptr<DispersionPolynomial>
    d(new DispersionPolynomial(disp_coeff, units::inv_cm, disp_var_1, ils_tab->band_name()));

  std::vector<boost::shared_ptr<Ils> > ils;
  ils.push_back(boost::shared_ptr<Ils>(new IlsGrating(d, ils_tab)));
  boost::shared_ptr<EmpiricalOrthogonalFunction> 
    eof(new EmpiricalOrthogonalFunction(1.0, *d, hf_eof, 0, 0, 1, ils_tab->band_name()));
  corr[0].push_back(eof);

  ils_tab.reset(new IlsTableLinear(hf_ils, 1, "WC-Band", "weak_co2"));

  disp_coeff = 5.74982835e+03, 1.99492886e-01;
  Array<double, 1> disp_var_2(3508);
  disp_var_2 = i1 + 1;

  d.reset(new DispersionPolynomial(disp_coeff, units::inv_cm, disp_var_2, ils_tab->band_name()));

  ils.push_back(boost::shared_ptr<Ils>(new IlsGrating(d, ils_tab)));

  ils_tab.reset(new IlsTableLinear(hf_ils, 2, "SC-Band", "strong_co2"));

  disp_coeff = 4.74980283e+03, 1.99492886e-01;
  Array<double, 1> disp_var_3(2005);
  disp_var_3 = i1 + 1;

  d.reset(new DispersionPolynomial(disp_coeff, units::inv_cm, disp_var_3, ils_tab->band_name()));
  ils.push_back(boost::shared_ptr<Ils>(new IlsGrating(d, ils_tab)));

  IlsInstrument inst(ils, corr);
  StateVector sv;
  sv.add_observer(inst);
  Array<double,1> x(5);
  x = 1.28695614e+04, 5.74982835e+03,4.74980283e+03,1.0,0;
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
  for(int i = 0; i < (int) plist.size(); ++i)
    rad_out_expect(i) += eof->eof().value(plist[i]);

  SpectralDomain sd(wn_in, units::inv_cm);
  SpectralRange sr(rad_hres_in, Unit("W / cm^2 / sr / cm^-1")); 
  BOOST_CHECK_MATRIX_CLOSE(inst.apply_instrument_model(Spectrum(sd, sr),
                                                       plist, 0).
                           spectral_range().data(), rad_out_expect);

  Array<double, 2> jac_rad_fake(rad_hres_in.rows(), 5);
  jac_rad_fake = 0;
  jac_rad_fake(Range::all(), 4) = rad_hres_in;
  ArrayAd<double, 1> rad_hres_in2(rad_hres_in, jac_rad_fake);
  SpectralRange sr2(rad_hres_in2, Unit("W / cm^2 / sr / cm^-1")); 
  BOOST_CHECK_MATRIX_CLOSE(inst.apply_instrument_model(Spectrum(sd, sr2),
                                                       plist, 0).
                           spectral_range().data(), rad_out_expect);
  Array<double, 2> jac = inst.apply_instrument_model(Spectrum(sd, sr2),
                                                     plist, 0).
    spectral_range().data_ad().jacobian();
  Array<double, 1> v0 = inst.apply_instrument_model(Spectrum(sd, sr2),
                                                    plist, 0).
    spectral_range().data();
  double epsilon = 1e-3;
  x(0) += epsilon;
  sv.update_state(x);
  Array<double, 1> v1 = inst.apply_instrument_model(Spectrum(sd, sr2),
                                                    plist, 0).
    spectral_range().data();
  Array<double, 2> jacd(jac.shape());
  jacd = 0;
  jacd(Range::all(), 0) = (v1 - v0) / epsilon;
  x(0) -= epsilon;
  x(3) += epsilon;
  sv.update_state(x);
  Array<double, 1> v2 = inst.apply_instrument_model(Spectrum(sd, sr2),
                                                    plist, 0).
    spectral_range().data();
  jacd(Range::all(), 3) = (v2 - v0) / epsilon;
  x(3) -= epsilon;
  sv.update_state(x);
  rad_hres_in2.value() *= 1 + epsilon;
  Array<double, 1> v3 = inst.apply_instrument_model(Spectrum(sd, sr2),
                                                    plist, 0).
    spectral_range().data();
  jacd(Range::all(), 4) = (v3 - v0) / epsilon;
  BOOST_CHECK_MATRIX_CLOSE(jac, jacd);
}
BOOST_AUTO_TEST_SUITE_END()


