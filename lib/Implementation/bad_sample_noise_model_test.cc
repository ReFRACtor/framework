#include "unit_test_support.h"
#include "bad_sample_noise_model.h"
#include "example_level_1b.h"
#include "fp_exception.h"
#include <iostream>

using namespace FullPhysics;
using namespace blitz;

class PassThroughNoiseModel : public NoiseModel {
public:
    PassThroughNoiseModel(const blitz::Array<double, 1> uncert_in) : uncert(uncert_in) {}
    blitz::Array<double, 1> uncertainty(int Spec_index, const blitz::Array<double, 1>& Radiance) const
        { return uncert; }
private:
    blitz::Array<double, 1> uncert;
};

BOOST_FIXTURE_TEST_SUITE(bad_sample_noise_model, GlobalFixture)
BOOST_AUTO_TEST_CASE(basic)
{
  boost::shared_ptr<HdfFile> h(new HdfFile(test_data_dir() + "in/common/l1b_example_data.h5"));
  ExampleLevel1b l1b = ExampleLevel1b(h, "2014090915251774");
  Array<double,1> uncertainty_calc = l1b.radiance(0).uncertainty();
  Array<bool, 2> bad_sample(3, uncertainty_calc.rows());
  bad_sample = false;
  bad_sample(0, 10) = true;
  boost::shared_ptr<NoiseModel> noise_model(new PassThroughNoiseModel(uncertainty_calc));
  BadSampleNoiseModel bnm(noise_model, bad_sample, 1e20);
  uncertainty_calc(10) = 1e20;
  BOOST_CHECK_MATRIX_CLOSE_TOL(bnm.uncertainty(0, l1b.radiance(0).data()), uncertainty_calc, 1e-3);
}

BOOST_AUTO_TEST_SUITE_END()

