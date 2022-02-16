#include "unit_test_support.h"

#include "atmosphere_fixture.h"
#include "stokes_coefficient_constant.h"
#include "altitude.h"

#include "pca_binning.h"
#include "pca_rt.h"

using namespace FullPhysics;
using namespace blitz;

void test_pca_rt(const boost::shared_ptr<AtmosphereStandard>& atm, const SpectralDomain& full_grid)
{
    bool debug_output = false;

    int num_streams = 8;
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

    // Sample the band away from the band continuum edge
    int beg_pca_range = 3000;
    int end_pca_range = 3999;
    if (full_grid.data().rows() <= end_pca_range) {
        Exception err;
        err << "PCA test range " << beg_pca_range << ", " << end_pca_range << " "
            << "exceeds side of full grid: " << full_grid.data().rows();
        throw err;
    }

    SpectralDomain test_grid( full_grid.data()(Range(beg_pca_range, end_pca_range)), full_grid.units() );

    // Use LIDORT for comparison
    LidortRt lidort_rt (atm, stokes_coefs, sza, zen, azm, pure_nadir,
                        num_streams, num_moments, false);
    Range r_all = Range::all();

    // Radiance alone
    blitz::Array<double, 2> pca_stokes = pca_rt.stokes(test_grid, 0);

    blitz::Array<double, 2> lidort_stokes = lidort_rt.stokes(test_grid, 0);

    BOOST_CHECK_MATRIX_CLOSE_TOL(lidort_stokes(r_all, 0), pca_stokes(r_all, 0), 2e-5);

    // Radiance and jacobian
    ArrayAd<double, 2> pca_stokes_and_jac = pca_rt.stokes_and_jacobian(test_grid, 0);

    ArrayAd<double, 2> lidort_stokes_and_jac = lidort_rt.stokes_and_jacobian(test_grid, 0);

    if(debug_output) {
        std::ofstream outfile_lid("pca_rt_test-lid_refl.out");
        outfile_lid << std::scientific
                    << std::setprecision(20)
                    << "# lidort refl" << std::endl
                    << lidort_stokes_and_jac.value()(r_all, 0) << std::endl;
        outfile_lid.close();
        std::ofstream outfile_pca("pca_rt_test-pca_refl.out");
        outfile_pca << std::scientific
                    << std::setprecision(20)
                    << "# pca refl" << std::endl
                    << pca_stokes_and_jac.value()(r_all, 0) << std::endl;
        outfile_pca.close();
    }

    BOOST_CHECK_MATRIX_CLOSE_TOL(lidort_stokes_and_jac.value()(r_all, 0), pca_stokes_and_jac.value()(r_all, 0), 2e-5);

    // Compare the RMS difference because the PCA jacobians are more of an approximation than the radiances
    // But that is okay since they are just used to drive the retrieval
    Array<double, 2> lid_jac = lidort_stokes_and_jac.jacobian()(r_all, 0, r_all);
    Array<double, 2> pca_jac = pca_stokes_and_jac.jacobian()(r_all, 0, r_all);

    if(debug_output) {
        std::ofstream outfile_lid("pca_rt_test-lid_jac.out");
        outfile_lid << std::scientific
                    << std::setprecision(20)
                    << "# lidort jac" << std::endl
                    << lid_jac << std::endl;
        outfile_lid.close();
        std::ofstream outfile_pca("pca_rt_test-pca_jac.out");
        outfile_pca << std::scientific
                    << std::setprecision(20)
                    << "# pca jac" << std::endl
                    << pca_jac << std::endl;
        outfile_pca.close();
    }

    firstIndex i1; secondIndex i2;
    Array<double, 2> jac_sqr_diff( sqr(lid_jac(i1, i2) - pca_jac(i1, i2)) );
    double rms_avg = 0;
    for(int jac_idx = 0; jac_idx < jac_sqr_diff.cols(); jac_idx++) {
        double rms = sqrt(mean(jac_sqr_diff(r_all, jac_idx)));
        BOOST_CHECK(rms < 0.2);
        if (debug_output) {
            std::cerr << "rms = " << rms << std::endl;
        }
        rms_avg += rms;
    }
    rms_avg /= jac_sqr_diff.cols();

    if (debug_output) {
        std::cerr << "rms average = " << rms_avg << std::endl;
    }

    BOOST_CHECK(rms_avg < 3e-3);
}

BOOST_FIXTURE_TEST_SUITE(pca_rt, AtmosphereFixture)

BOOST_AUTO_TEST_CASE(with_aerosol)
{
  test_pca_rt(atm, highres_grid(0));
}

BOOST_AUTO_TEST_CASE(no_aerosol)
{ 
  boost::shared_ptr<AtmosphereStandard> src_atm(boost::dynamic_pointer_cast<AtmosphereStandard>(atm));

  std::vector<boost::shared_ptr<Altitude> > alt_clone;
    BOOST_FOREACH(const boost::shared_ptr<Altitude>& a, src_atm->altitude_ptr())
        alt_clone.push_back(a->clone());

  boost::shared_ptr<AtmosphereStandard>
    atm_no_aerosol(new AtmosphereStandard(src_atm->absorber_ptr()->clone(),
                                          src_atm->pressure_ptr()->clone(),
                                          src_atm->temperature_ptr()->clone(),
                                          src_atm->rayleigh_ptr()->clone(),
                                          NULL,
                                          src_atm->relative_humidity_ptr()->clone(),
                                          src_atm->ground()->clone(),
                                          alt_clone,
                                          src_atm->constant_ptr()));

  test_pca_rt(atm_no_aerosol, highres_grid(0));
 
}

BOOST_AUTO_TEST_SUITE_END()
