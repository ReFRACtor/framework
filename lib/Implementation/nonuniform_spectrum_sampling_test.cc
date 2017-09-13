#include "nonuniform_spectrum_sampling.h"
#include "uniform_spectrum_sampling.h"
#include "unit_test_support.h"
#include "fp_exception.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(nonuniform_spectrum_sampling, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic_test)
{
  IfstreamCs grid_in1(test_data_dir() + "in/nonunif_rt_grid/nonunif_rt_grid_abo2");
  IfstreamCs grid_in2(test_data_dir() + "in/nonunif_rt_grid/nonunif_rt_grid_wco2");
  IfstreamCs grid_in3(test_data_dir() + "in/nonunif_rt_grid/nonunif_rt_grid_sco2");

  IfstreamCs grid_in3_sorted(test_data_dir() + "in/nonunif_rt_grid/nonunif_rt_grid_sco2_sorted");

  Array<double, 1> grid_data1, grid_data2, grid_data3, grid_data3_sorted;
  grid_in1 >> grid_data1;
  grid_in2 >> grid_data2;
  grid_in3 >> grid_data3;
  grid_in3_sorted >> grid_data3_sorted;

  SpectralDomain grid_dom1(grid_data1, units::inv_cm);
  SpectralDomain grid_dom2(grid_data2, units::inv_cm);
  SpectralDomain grid_dom3(grid_data3, units::inv_cm);
  SpectralDomain grid_dom3_sorted(grid_data3_sorted, units::inv_cm);

  boost::shared_ptr<SpectrumSampling> 
    interpolated_spec(new UniformSpectrumSampling(31.8,34.5,0.05,
                                                  1.5, 10.5, 0.05,
                                                  21.7, 26.3, 0.05));
  boost::shared_ptr<SpectrumSampling> 
    interpolated_spec1(new UniformSpectrumSampling(31.8,34.5,0.05));
  boost::shared_ptr<SpectrumSampling> 
    interpolated_spec2(new UniformSpectrumSampling(1.5, 10.5, 0.05));
  boost::shared_ptr<SpectrumSampling> 
    interpolated_spec3(new UniformSpectrumSampling(21.7, 26.3, 0.05));
  NonuniformSpectrumSampling nonunif_rt_grid(grid_dom1, grid_dom2, grid_dom3, interpolated_spec);
  NonuniformSpectrumSampling nonunif_rt_grid_1(grid_dom1, interpolated_spec1);
  NonuniformSpectrumSampling nonunif_rt_grid_2(grid_dom2, interpolated_spec2);
  NonuniformSpectrumSampling nonunif_rt_grid_3_srtd(grid_dom3_sorted, interpolated_spec3);

  if(false) {                        // Can turn on debugging messages if desired.
    std::cout << nonunif_rt_grid << "\n"
              << nonunif_rt_grid_1 << "\n"
              << nonunif_rt_grid_2 << "\n"
              << nonunif_rt_grid_3_srtd << "\n";
  }

  SpectralDomain dummy;                // Not used by UniformSpectrumSampling
  DoubleWithUnit dummy2;
  BOOST_CHECK_EQUAL(nonunif_rt_grid.
                    spectral_domain(0, dummy, dummy2).wavenumber().size(), 7);
  BOOST_CHECK_EQUAL(nonunif_rt_grid.
                    spectral_domain(1, dummy, dummy2).wavenumber().size(), 10);
  BOOST_CHECK_EQUAL(nonunif_rt_grid.
                    spectral_domain(2, dummy, dummy2).wavenumber().size(), 9);

  BOOST_CHECK_EQUAL(nonunif_rt_grid_1.
                    spectral_domain(0, dummy, dummy2).wavenumber().size(), 7);
  BOOST_CHECK_EQUAL(nonunif_rt_grid_2.
                    spectral_domain(0, dummy, dummy2).wavenumber().size(), 10);
  BOOST_CHECK_EQUAL(nonunif_rt_grid_3_srtd.
                    spectral_domain(0, dummy, dummy2).wavenumber().size(), 9);

  BOOST_CHECK_EQUAL(sum(abs(nonunif_rt_grid.
                    spectral_domain(2, dummy, dummy2).wavenumber()-
                    nonunif_rt_grid_3_srtd.
                    spectral_domain(0, dummy, dummy2).wavenumber())), 0.0);
}

BOOST_AUTO_TEST_SUITE_END()
