#include "instrument_measurement_level_1b.h"

using namespace FullPhysics;
using namespace blitz;

InstrumentMeasurementLevel1b::InstrumentMeasurementLevel1b(const boost::shared_ptr<Level1b>& level_1b, const boost::shared_ptr<ForwardModelSpectralGrid>& spectral_grids)
    : l1b(level_1b), grids(spectral_grids)
{
}

int InstrumentMeasurementLevel1b::num_channels() const
{
    return grids->number_spectrometer();
}

const SpectralDomain InstrumentMeasurementLevel1b::spectral_grid(int channel_index) const
{
    return grids->low_resolution_grid(channel_index);
}

Spectrum InstrumentMeasurementLevel1b::radiance(int channel_index) const
{
    range_check(channel_index, 0, num_channels());
    Spectrum full(spectral_grid(channel_index), l1b->radiance(channel_index));
    const std::vector<int>& plist = grids->pixel_list(channel_index);
    Array<double, 1> res_d((int) plist.size());
    ArrayAd<double, 1> res_r((int) plist.size(), full.spectral_range().data_ad().number_variable());
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
                    SpectralRange(res_r, full.spectral_range().units(),
                                  uncer));
}
