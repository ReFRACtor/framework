#include "cloud_3d_effect.h"
#include "unit_test_support.h"
#include "global_fixture.h"

#include "forward_model_spectral_grid.h"

using namespace FullPhysics;
using namespace blitz;

class Cloud3dFixture : public GlobalFixture {
public:
    Cloud3dFixture()
    {
        initial_offset = 0.2;;
        initial_slope = 0.0002;
        band_name = "TEST";

        cloud_3d.reset(new Cloud3dEffect(initial_offset, initial_slope, band_name));
    }
    virtual ~Cloud3dFixture() = default;

    double initial_offset;
    double initial_slope;
    std::string band_name;
 
    boost::shared_ptr<Cloud3dEffect> cloud_3d;
};

BOOST_FIXTURE_TEST_SUITE(cloud_3d_effect, Cloud3dFixture)

BOOST_AUTO_TEST_CASE(init)
{
    BOOST_CHECK_CLOSE(initial_offset, cloud_3d->offset().value(), 1e-10);
    BOOST_CHECK_CLOSE(initial_slope, cloud_3d->slope().value(), 1e-10);
}

BOOST_AUTO_TEST_CASE(apply)
{
    Array<double, 1> inp_vals(10);
    firstIndex i1;
    inp_vals = i1;

    Array<double, 1> expt_vals(inp_vals.shape());
    expt_vals = inp_vals(i1) * (1 + initial_offset + initial_slope * inp_vals(i1));

    Spectrum spec(SpectralDomain(inp_vals, units::inv_cm), SpectralRange(inp_vals, units::dimensionless));

    ForwardModelSpectralGrid ignored;
    cloud_3d->apply_effect(spec, ignored);

    BOOST_CHECK_MATRIX_CLOSE(expt_vals, spec.spectral_range().data());
}


BOOST_AUTO_TEST_SUITE_END()
