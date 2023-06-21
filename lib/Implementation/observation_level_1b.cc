#include "observation_level_1b.h"
#include "level_1b_sample_coefficient.h"

#include "fp_serialize_support.h"
using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void ObservationLevel1b::serialize(Archive& ar, const unsigned int UNUSED(version))
{
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Observation)
       & FP_NVP(l1b) & FP_NVP(inst) & FP_NVP(grids);
}

FP_IMPLEMENT(ObservationLevel1b);
#endif

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(ObservationLevel1b, Observation)
.def(luabind::constructor<const boost::shared_ptr<Level1b>&, const boost::shared_ptr<Instrument> &, const boost::shared_ptr<ForwardModelSpectralGrid>&>())
.def(luabind::constructor<const boost::shared_ptr<Level1bSampleCoefficient>&, const boost::shared_ptr<Instrument> &, const boost::shared_ptr<ForwardModelSpectralGrid>&>())
REGISTER_LUA_END()
#endif

ObservationLevel1b::ObservationLevel1b(const boost::shared_ptr<Level1b>& level_1b,
        const boost::shared_ptr<Instrument> &instrument,
        const boost::shared_ptr<ForwardModelSpectralGrid>& spectral_grids)
    : l1b(level_1b), inst(instrument), grids(spectral_grids)
{
}

int ObservationLevel1b::num_channels() const
{
    return grids->number_spectrometer();
}

SpectralDomain ObservationLevel1b::spectral_domain(int sensor_index) const
{
    return grids->low_resolution_grid(sensor_index);
}

Spectrum ObservationLevel1b::radiance(int sensor_index, bool skip_jacobian) const
{
    range_check(sensor_index, 0, num_channels());

    Spectrum full(inst->pixel_spectral_domain(sensor_index), l1b->radiance(sensor_index));

    const std::vector<int>& plist = grids->pixel_list(sensor_index);
    Array<double, 1> res_d((int) plist.size());

    int num_jac;
    if (skip_jacobian) {
        num_jac = 0;
    } else {
        num_jac = full.spectral_range().data_ad().number_variable();
    }

    ArrayAd<double, 1> res_r((int) plist.size(), num_jac);
    Array<double, 1> uncer;

    if(full.spectral_range().uncertainty().rows() > 0) {
        uncer.resize(res_r.rows());
    }

    for(int i = 0; i < res_d.rows(); ++i) {
        res_d(i) = full.spectral_domain().data()(plist[i]);
        AutoDerivative<double> t = full.spectral_range().data_ad()(plist[i]);
        res_r(i) = t;

        if(uncer.rows() > 0) {
            uncer(i) = full.spectral_range().uncertainty()(plist[i]);
        }
    }

    return Spectrum(SpectralDomain(res_d, full.spectral_domain().units()),
                    SpectralRange(res_r, full.spectral_range().units(), uncer));
}
