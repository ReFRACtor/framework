#include "sample_grid_spectral_domain.h"
#include "fp_serialize_support.h"
#include <boost/lexical_cast.hpp>
using namespace FullPhysics;
using namespace blitz;

#ifdef FP_HAVE_BOOST_SERIALIZATION

template<class Archive>
void SampleGridSpectralDomain::serialize(Archive & ar,
                                  const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(SampleGridImpBase)
    & FP_NVP(spec_domain) & FP_NVP_(band_name);
}

FP_IMPLEMENT(SampleGridSpectralDomain);
#endif

#ifdef HAVE_LUA
#include "register_lua.h"

REGISTER_LUA_DERIVED_CLASS(SampleGridSpectralDomain, SampleGrid)
.def(luabind::constructor<const SpectralDomain&, const std::string&>())
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Constructor.
//-----------------------------------------------------------------------

SampleGridSpectralDomain::SampleGridSpectralDomain(const SpectralDomain& Spec_domain, const std::string& Band_name)
    : band_name_(Band_name)
{
    // Create index array for spectral domain if not already present
    if (Spec_domain.sample_index().rows() == 0) {
        blitz::Array<int, 1> spectral_index(Spec_domain.data().rows());
        for(int i = 0; i < spectral_index.rows(); ++i) {
            spectral_index(i) = i + 1;
        }
        spec_domain = SpectralDomain(Spec_domain.data(), spectral_index, Spec_domain.units());
    }
    else {
        spec_domain = Spec_domain;
    }
}

// See base class for description.
SpectralDomain
SampleGridSpectralDomain::sample_grid() const
{
    return spec_domain;
}

boost::shared_ptr<SampleGrid> SampleGridSpectralDomain::clone() const
{
    return boost::shared_ptr<SampleGrid>(new SampleGridSpectralDomain(spec_domain, band_name_));
}

void SampleGridSpectralDomain::print(std::ostream& Os) const
{
    Os << "SampleGridSpectralDomain for band " << band_name_ << "\n"
       << "  Spectral Domain:  " << spec_domain.data() << "\n";
}
