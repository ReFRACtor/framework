#include "sample_grid_spectral_domain.h"
#include <boost/lexical_cast.hpp>
using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"

REGISTER_LUA_DERIVED_CLASS(SampleGridSpectralDomain, SampleGrid)
.def(luabind::constructor<const SpectralDomain&, const std::string&,
                          int, bool>())
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Constructor.
//-----------------------------------------------------------------------

SampleGridSpectralDomain::SampleGridSpectralDomain
(const SpectralDomain& Spec_domain,
 const std::string& Band_name,
 int Number_pixel,
 bool Is_one_based)
: spec_domain(Spec_domain),
  is_one_based(Is_one_based),
  band_name_(Band_name),
  index_array(Number_pixel),
  spectral_index(Number_pixel)
{
  initialize();
}



// Initialize class internals
void SampleGridSpectralDomain::initialize() {
  for(int i = 0; i < index_array.rows(); ++i) {
    index_array(i) = i + (is_one_based ? 1 : 0);
    spectral_index(i) = i + 1;
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
  return boost::shared_ptr<SampleGrid>
    (new SampleGridSpectralDomain(spec_domain, band_name_,
                                  index_array.rows(),
			                      is_one_based));
}

void SampleGridSpectralDomain::print(std::ostream& Os) const
{
  Os << "SampleGridSpectralDomain for band " << band_name_ << "\n"
     << "  1 based:  " << (is_one_based ? "True" : "False") << "\n"
     << "  Spectral Domain:  " << spec_domain.data() << "\n";
}
