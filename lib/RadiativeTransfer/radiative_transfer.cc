#include "radiative_transfer.h"
#include "fp_serialize_support.h"
#include "logger.h"
using namespace FullPhysics;

AccumulatedTimer RadiativeTransfer::timer("RT");

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void RadiativeTransfer::serialize(Archive & ar,
			const unsigned int UNUSED(version))
{
  FP_GENERIC_BASE(RadiativeTransfer);

  // Dummy placeholder, just so we can have derived classes call
  // serialization of this. We use to have derived classes "know"
  // that the base class doesn't have anything. But seems better to
  // *always* have base classes do something, so we can add stuff in
  // the future w/o breaking a bunch of code.
  std::string p = "empty";
  ar & FP_NVP2("placeholder", p);
}

FP_IMPLEMENT(RadiativeTransfer);
#endif

using namespace FullPhysics;
#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_CLASS(RadiativeTransfer)
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Helper routine, creates a progress meter. This will return 0 if we
/// aren't logging, or if we don't have enough points to bother with.
//-----------------------------------------------------------------------
boost::shared_ptr<boost::progress_display> 
RadiativeTransfer::progress_display(const blitz::Array<double, 1>& wn, boost::optional<std::string> message) const
{
  boost::shared_ptr<boost::progress_display> res;
  if(wn.size() > 100 && Logger::stream()) {
    if (!message) {
        Logger::info() << "RT Progress\n";
    } else {
        Logger::info() << *message << "\n";
    }
    res.reset(new boost::progress_display(wn.size(), *Logger::stream()));
  }
  return res;
}
