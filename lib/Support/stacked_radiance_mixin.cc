#include "stacked_radiance_mixin.h"
#include "fp_serialize_support.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void StackedRadianceMixin::serialize(Archive & ar,
				     const unsigned int UNUSED(version))
{
  FP_GENERIC_BASE(StackedRadianceMixin);

  // Dummy placeholder, just so we can have derived classes call
  // serialization of this. We use to have derived classes "know"
  // that the base class doesn't have anything. But seems better to
  // *always* have base classes do something, so we can add stuff in
  // the future w/o breaking a bunch of code.
  std::string p = "empty";
  ar & FP_NVP2("placeholder", p);
}

FP_IMPLEMENT(StackedRadianceMixin);
#endif

boost::optional<Range> StackedRadianceMixin::stacked_pixel_range(int channel_index) const
{
    range_check(channel_index, 0, num_channels());
    int sind = 0;

    for(int i = 0; i < channel_index; ++i) {
        sind += spectral_domain(i).data().rows();
    }

    int nrow = spectral_domain(channel_index).data().rows();

    if(nrow > 0) {
        return boost::optional<Range>(Range(sind, sind + nrow - 1));
    } else {
        return boost::optional<Range>();
    }
}


Spectrum StackedRadianceMixin::radiance_all(bool skip_jacobian) const
{
    std::vector<Spectrum> sall;
    std::vector<Range> prall;

    for(int i = 0; i < num_channels(); ++i) {
        boost::optional<Range> pr = stacked_pixel_range(i);

        if(pr) {
            sall.push_back(radiance(i, skip_jacobian));
            prall.push_back(*pr);
        }
    }

    if(sall.size() == 0) {
        throw Exception("Measured radiance Spectrum empty, pixel ranges must be empty");
    }

    int nrow = 0;
    int nvar = 0;
    bool have_uncertainty = true;
    BOOST_FOREACH(const Spectrum & s, sall) {
        nrow += s.spectral_domain().data().rows();
        nvar = std::max(nvar, s.spectral_range().data_ad().number_variable());

        if(s.spectral_range().uncertainty().rows() == 0) {
            have_uncertainty = false;
        }
    }
    Unit ud, ur;

    if(sall.size() > 0) {
        ud = sall[0].spectral_domain().units();
        ur = sall[0].spectral_range().units();
    }

    ArrayAd<double, 1> sr(nrow, nvar);
    Array<double, 1> sd(nrow);
    Array<double, 1> uncer;

    if(have_uncertainty) {
        uncer.resize(nrow);
    }

    for(int i = 0; i < (int) sall.size(); ++i) {
        sd(prall[i]) = sall[i].spectral_domain().data() *
                       FullPhysics::conversion(sall[i].spectral_domain().units(), ud);
        sr.value()(prall[i]) = sall[i].spectral_range().data() *
                               FullPhysics::conversion(sall[i].spectral_range().units(), ur);

        if(have_uncertainty)
            uncer(prall[i]) = sall[i].spectral_range().uncertainty() *
                              FullPhysics::conversion(sall[i].spectral_range().units(), ur);

        if(nvar > 0)
            sr.jacobian()(prall[i], Range::all()) =
                sall[i].spectral_range().data_ad().jacobian() *
                FullPhysics::conversion(sall[i].spectral_range().units(), ur);
    }

    return Spectrum(SpectralDomain(sd, ud), SpectralRange(sr, ur, uncer));
}
