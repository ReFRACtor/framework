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

    blitz::Array<double, 2> pca_stokes = pca_rt.stokes(test_grid, 0);

    LidortRt lidort_rt (atm, stokes_coefs, sza, zen, azm, pure_nadir,
                        num_streams, num_moments, false);

    blitz::Array<double, 2> lidort_stokes = lidort_rt.stokes(test_grid, 0);

    Range r_all = Range::all();
    BOOST_CHECK_MATRIX_CLOSE_TOL(lidort_stokes(r_all, 0), pca_stokes(r_all, 0), 1e-6);
   
}

BOOST_AUTO_TEST_SUITE_END()
