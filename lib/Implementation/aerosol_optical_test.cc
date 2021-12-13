#include "aerosol_optical.h"
#include "pressure.h"
#include "unit_test_support.h"
#include "atmosphere_fixture.h"
#include "aerosol_property_hdf.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(aerosol_optical, AtmosphereFixture)

BOOST_AUTO_TEST_CASE(config_file)
{
  boost::shared_ptr<AerosolOptical> a =
    boost::dynamic_pointer_cast<AerosolOptical>(config_aerosol);
  Array<double, 1> od_expect(4);
  od_expect = 0.0124999579334943, 0.0124999579334943, 0.012499957933494298, 0.0124999579334943;
  for(int aer_idx = 0; aer_idx < od_expect.rows(); aer_idx++)
    BOOST_CHECK_CLOSE(a->aerosol_optical_depth(aer_idx), od_expect(aer_idx), 1e-4);
  double total_expect = sum(od_expect);
  BOOST_CHECK_CLOSE(a->aerosol_optical_depth_total(), total_expect, 1e-4);
}

BOOST_AUTO_TEST_CASE(layer_parameters)
{
  firstIndex i1; secondIndex i2;
  // Expected values were gotten by running the old Fortran code and
  // extracting out the answer from that.
  boost::shared_ptr<AerosolOptical> a =
    boost::dynamic_pointer_cast<AerosolOptical>(config_aerosol);

  IfstreamCs expt_input(test_data_dir() + "expected/aerosol_optical/layer_parameters");

  Array<double, 1> ext_od_expect;
  Array<double, 1> sca_od_expect;

  expt_input >> ext_od_expect >> sca_od_expect;

  Array<double, 1> ext_od_calc(sum(a->extinction_optical_depth_each_layer(12930.30).value(), i2));
  Array<double, 1> sca_od_calc(sum(a->scattering_optical_depth_each_layer(12930.30).value(), i2));

  BOOST_CHECK_MATRIX_CLOSE(ext_od_calc, ext_od_expect);
  BOOST_CHECK_MATRIX_CLOSE(sca_od_calc, sca_od_expect);
  
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(aerosol_jac, LuaConfigurationFixture)

BOOST_AUTO_TEST_CASE(optical_depth_jac)
{
  boost::shared_ptr<AerosolOptical> a =
    boost::dynamic_pointer_cast<AerosolOptical>(config_aerosol);
  StateVector& sv = *config_state_vector;
  Array<double, 1> sv0(sv.state().copy());
  double wn = 4820.0;
  ArrayAd<double, 2> od = a->extinction_optical_depth_each_layer(wn);
  Array<double, 2> od0(od.shape());
  od0 = od.value();
  Array<double, 3> jac = od.jacobian().copy();
  for(int i = 0; i < sv.state().rows(); ++i) {
    Array<double, 1> svn(sv0.copy());
    svn(i) += epsilon(i);
    sv.update_state(svn);
    Array<double, 2> jacfd(od0.shape());
    jacfd = (a->extinction_optical_depth_each_layer(wn).value() - od0) / epsilon(i);
    if(false) {                        // Can turn this off to dump values,
                                // if needed for debugging
      double diff = max(abs(jac(Range::all(), Range::all(), i) - jacfd));
      if(diff > 0)
        std::cerr << i << ": " << jac(Range::all(), Range::all(), i) << "\n"
                  << jacfd << "\n"
                  << diff << "\n";
    }
    BOOST_CHECK_MATRIX_CLOSE_TOL(jac(Range::all(), Range::all(), i), jacfd, 1e-6);
  }
}

BOOST_AUTO_TEST_CASE(sca_od_jac)
{
  boost::shared_ptr<AerosolOptical> a =
    boost::dynamic_pointer_cast<AerosolOptical>(config_aerosol);
  StateVector& sv = *config_state_vector;
  Array<double, 1> sv0(sv.state().copy());
  double wn = 4820.0;
  ArrayAd<double, 2> sca_od = a->scattering_optical_depth_each_layer(wn);
  Array<double, 2> sca_od0(sca_od.shape());
  sca_od0 = sca_od.value();
  Array<double, 3> jac = sca_od.jacobian().copy();
  for(int i = 0; i < sv.state().rows(); ++i) {
    Array<double, 1> svn(sv0.copy());
    svn(i) += epsilon(i);
    sv.update_state(svn);
    Array<double, 2> jacfd(sca_od0.shape());
    jacfd = (a->scattering_optical_depth_each_layer(wn).value() - sca_od0) / epsilon(i);
    if(false) {                        // Can turn this off to dump values,
                                // if needed for debugging
      double diff = max(abs(jac(Range::all(), Range::all(), i) - jacfd));
      if(diff > 0)
        std::cerr << i << ": " << jac(Range::all(), Range::all(), i) << "\n"
                  << jacfd << "\n"
                  << diff << "\n";
    }
    BOOST_CHECK_MATRIX_CLOSE_TOL(jac(Range::all(), Range::all(), i), jacfd, 1e-6);
  }
}

BOOST_AUTO_TEST_SUITE_END()


