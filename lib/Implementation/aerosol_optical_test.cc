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

  Array<double, 1> od_expect;
  Array<double, 1> ssa_expect;

  expt_input >> od_expect >> ssa_expect;

  Array<double, 1> od_calc(sum(a->optical_depth_each_layer(12930.30).value(), i2));
  Array<double, 1> ssa_calc(a->ssa_each_layer(12930.30).value());

  BOOST_CHECK_MATRIX_CLOSE(od_calc, od_expect);
  BOOST_CHECK_MATRIX_CLOSE(ssa_calc, ssa_expect);
  
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(aerosol_jac, ConfigurationFixture)

BOOST_AUTO_TEST_CASE(optical_depth_jac)
{
  boost::shared_ptr<AerosolOptical> a =
    boost::dynamic_pointer_cast<AerosolOptical>(config_aerosol);
  StateVector& sv = *config_state_vector;
  Array<double, 1> sv0(sv.state().copy());
  double wn = 4820.0;
  ArrayAd<double, 2> od = a->optical_depth_each_layer(wn);
  Array<double, 2> od0(od.shape());
  od0 = od.value();
  Array<double, 3> jac = od.jacobian().copy();
  for(int i = 0; i < sv.state().rows(); ++i) {
    Array<double, 1> svn(sv0.copy());
    svn(i) += epsilon(i);
    sv.update_state(svn);
    Array<double, 2> jacfd(od0.shape());
    jacfd = (a->optical_depth_each_layer(wn).value() - od0) / epsilon(i);
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

BOOST_AUTO_TEST_CASE(ssa_jac)
{
  boost::shared_ptr<AerosolOptical> a =
    boost::dynamic_pointer_cast<AerosolOptical>(config_aerosol);
  StateVector& sv = *config_state_vector;
  Array<double, 1> sv0(sv.state().copy());
  double wn = 4820.0;
  ArrayAd<double, 1> ssa = a->ssa_each_layer(wn);
  Array<double, 1> ssa0(ssa.shape());
  ssa0 = ssa.value();
  Array<double, 2> jac = ssa.jacobian().copy();
  for(int i = 0; i < sv.state().rows(); ++i) {
    Array<double, 1> svn(sv0.copy());
    svn(i) += epsilon(i);
    sv.update_state(svn);
    Array<double, 1> jacfd(ssa0.shape());
    jacfd = (a->ssa_each_layer(wn).value() - ssa0) / epsilon(i);
    if(false) {                        // Can turn this off to dump values,
                                // if needed for debugging
      double diff = max(abs(jac(Range::all(), i) - jacfd));
      if(diff > 0)
        std::cerr << i << ": " << jac(Range::all(), i) << "\n"
                  << jacfd << "\n"
                  << diff << "\n";
    }
    BOOST_CHECK_MATRIX_CLOSE_TOL(jac(Range::all(), i), jacfd, 1e-6);
  }
}

BOOST_AUTO_TEST_SUITE_END()


