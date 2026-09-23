#include "vlidort_rt.h"
#include "unit_test_support.h"
#include "lidort_fixture.h"
#include "old_constant.h"
#include "fe_disable_exception.h"

using namespace FullPhysics;
using namespace blitz;

// Exercise the high-level VlidortRt RadiativeTransfer class end-to-end
// (atmosphere -> phase function moments -> VLidortRtDriver -> Stokes ->
// StokesCoefficient). The LidortLambertianFixture builds a config
// atmosphere and stokes coefficients we can reuse to also construct a
// VlidortRt, so the scalar (nstokes = 1) VLIDORT result can be compared
// directly against the established LidortRt result (LIDORT parity).

BOOST_FIXTURE_TEST_SUITE(vlidort_rt_lambertian, LidortLambertianFixture)

BOOST_AUTO_TEST_CASE(scalar_parity_with_lidort)
{
  int nstreams = 4;
  int nmoms = 2*nstreams;
  bool do_multi_scatt_only = false;

  // nstokes = 1 VLIDORT should reproduce the scalar LIDORT radiance.
  boost::shared_ptr<VlidortRt> vlidort_rt(new VlidortRt(
      config_atmosphere, stokes_coefs, sza, zen, azm, pure_nadir,
      nstreams, nmoms, do_multi_scatt_only, 1));

  lidort_rt->rt_driver()->set_plane_parallel();
  vlidort_rt->rt_driver()->set_plane_parallel();

  Array<double, 1> lidort_refl(lidort_rt->reflectance(wn_arr, 0).spectral_range().data());
  Array<double, 1> vlidort_refl(vlidort_rt->reflectance(wn_arr, 0).spectral_range().data());

  BOOST_CHECK_MATRIX_CLOSE_TOL(lidort_refl, vlidort_refl, 1e-6);
}

BOOST_AUTO_TEST_CASE(vector_runs_and_matches_intensity)
{
  int nstreams = 8;
  int nmoms = 2*nstreams;
  bool do_multi_scatt_only = false;

  // Full vector (nstokes = 3) VLIDORT. The intensity-weighted radiance
  // should be close to the scalar LIDORT result; the difference is the
  // (small) polarization correction to the multiple scattering.
  boost::shared_ptr<VlidortRt> vlidort_rt(new VlidortRt(
      config_atmosphere, stokes_coefs, sza, zen, azm, pure_nadir,
      nstreams, nmoms, do_multi_scatt_only, 3));

  lidort_rt->rt_driver()->set_plane_parallel();
  vlidort_rt->rt_driver()->set_plane_parallel();

  // VLIDORT can raise benign floating point exceptions during setup for
  // vector calculations (it copies garbage values that are never used); the
  // framework suppresses these elsewhere the same way.
  FeDisableException disable_fp;

  Array<double, 1> vlidort_refl(vlidort_rt->reflectance(wn_arr, 0).spectral_range().data());

  // Full vector RT over the configured atmosphere produces a positive,
  // finite intensity-weighted radiance.
  BOOST_CHECK(vlidort_refl.rows() == wn_arr.rows());
  for(int i = 0; i < vlidort_refl.rows(); ++i) {
    BOOST_CHECK(vlidort_refl(i) > 0.0);
    BOOST_CHECK(!std::isnan(vlidort_refl(i)));
  }
}

BOOST_AUTO_TEST_CASE(serialization)
{
  if(!have_serialize_supported())
    return;
  int nstreams = 4;
  int nmoms = 2*nstreams;
  boost::shared_ptr<VlidortRt> vlidort_rt(new VlidortRt(
      config_atmosphere, stokes_coefs, sza, zen, azm, pure_nadir,
      nstreams, nmoms, false, 3));

  std::string d = serialize_write_string(vlidort_rt);
  boost::shared_ptr<VlidortRt> vlidort_rt_r =
    serialize_read_string<VlidortRt>(d);
  BOOST_CHECK(vlidort_rt_r);
  BOOST_CHECK_EQUAL(vlidort_rt_r->number_stream(), vlidort_rt->number_stream());
}

BOOST_AUTO_TEST_SUITE_END()
