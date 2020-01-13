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

  SpectralDomain sd = highres_grid(0);

  blitz::Array<double, 2> stokes = pca_rt.stokes(sd, 0);
 
}

BOOST_AUTO_TEST_SUITE_END()
