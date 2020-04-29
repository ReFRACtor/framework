#include "unit_test_support.h"

#include "atmosphere_fixture.h"
#include "stokes_coefficient_constant.h"
#include "pca_binning.h"

#include "pca_rt.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(pca_rt, AtmosphereFixture)

BOOST_AUTO_TEST_CASE(basic)
{
    int num_streams = 4;
    int num_moments = 2 * num_streams;
    bool pure_nadir = false;

    std::string primary_absorber = "CO2";

    PCABinning::Method bin_method = PCABinning::Method::UVVSWIR_V4;
    int num_bins = 11;
    int num_eofs = 4;

    blitz::Array<double, 2> stokes_coef_v(3, 3);
    stokes_coef_v = 
        1,0,0,
        1,0,0,
        1,0,0;
    boost::shared_ptr<StokesCoefficientConstant> stokes_coefs(new StokesCoefficientConstant(stokes_coef_v));

    int n_channels = 3;
    blitz::Array<double, 1> sza(n_channels), zen(n_channels), azm(n_channels);
    sza = 0.001;
    zen = 0.001;
    azm = 0.0;

    PCARt pca_rt (atm, primary_absorber,
                  bin_method, num_bins, num_eofs,
                  stokes_coefs,
                  sza, zen, azm,
                  num_streams, num_moments);

    SpectralDomain full_grid = highres_grid(0);
    SpectralDomain test_grid( full_grid.data()(Range(0, 999)), full_grid.units() );

    // Use LIDORT for comparison
    LidortRt lidort_rt (atm, stokes_coefs, sza, zen, azm, pure_nadir,
                        num_streams, num_moments, false);
    Range r_all = Range::all();

    // Radiance alone
    blitz::Array<double, 2> pca_stokes = pca_rt.stokes(test_grid, 0);

    blitz::Array<double, 2> lidort_stokes = lidort_rt.stokes(test_grid, 0);

    BOOST_CHECK_MATRIX_CLOSE_TOL(lidort_stokes(r_all, 0), pca_stokes(r_all, 0), 1e-6);

    // Radiance and jacobian
    ArrayAd<double, 2> pca_stokes_and_jac = pca_rt.stokes_and_jacobian(test_grid, 0);

    ArrayAd<double, 2> lidort_stokes_and_jac = lidort_rt.stokes_and_jacobian(test_grid, 0);

    BOOST_CHECK_MATRIX_CLOSE_TOL(lidort_stokes_and_jac.value()(r_all, 0), pca_stokes_and_jac.value()(r_all, 0), 1e-6);

    // Compare the RMS difference because the PCA jacobians are more of an approximation than the radiances
    // But that is okay since they are just used to drive the retrieval
    Array<double, 2> lid_jac = lidort_stokes_and_jac.jacobian()(r_all, 0, r_all);
    Array<double, 2> pca_jac = pca_stokes_and_jac.jacobian()(r_all, 0, r_all);

    firstIndex i1; secondIndex i2;
    Array<double, 2> jac_sqr_diff( sqr(lid_jac(i1, i2) - pca_jac(i1, i2)) );
    double rms_avg = 0;
    for(int jac_idx = 0; jac_idx < jac_sqr_diff.cols(); jac_idx++) {
        double rms = sqrt(mean(jac_sqr_diff(r_all, jac_idx)));
        BOOST_CHECK(rms < 0.7);
        rms_avg += rms;
    }
    rms_avg /= jac_sqr_diff.cols();

    BOOST_CHECK(rms_avg < 0.02);
}

BOOST_AUTO_TEST_SUITE_END()
