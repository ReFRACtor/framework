#include "observation_level_1b.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(ObservationLevel1b, Observation)
.def(luabind::constructor<const boost::shared_ptr<Level1bSampleCoefficient>&, const boost::shared_ptr<Instrument> &, const boost::shared_ptr<ForwardModelSpectralGrid>&>())
REGISTER_LUA_END()
#endif

ObservationLevel1b::ObservationLevel1b(const boost::shared_ptr<Level1bSampleCoefficient>& level_1b,
        const boost::shared_ptr<Instrument> &instrument,
        const boost::shared_ptr<ForwardModelSpectralGrid>& spectral_grids)
    : l1b(level_1b), inst(instrument), grids(spectral_grids)
{
}

int ObservationLevel1b::num_channels() const
{
    return grids->number_spectrometer();
}

const SpectralDomain ObservationLevel1b::spectral_domain(int channel_index) const
{
    return grids->low_resolution_grid(channel_index);
}

Spectrum ObservationLevel1b::radiance(int channel_index, bool skip_jacobian) const
{
    range_check(channel_index, 0, num_channels());

    Spectrum full(inst->pixel_spectral_domain(channel_index), l1b->radiance(channel_index));

    const std::vector<int>& plist = grids->pixel_list(channel_index);
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
