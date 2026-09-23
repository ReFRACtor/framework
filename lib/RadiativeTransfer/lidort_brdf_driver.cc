#include "lidort_brdf_driver.h"

#include "fp_serialize_support.h"
#include "fp_exception.h"
#include "linear_algebra.h"
#include "ostream_pad.h"
#include "fe_disable_exception.h"
#include "fp_exception.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef FP_HAVE_BOOST_SERIALIZATION

template<class Archive>
void LidortBrdfDriver::serialize(Archive & ar, const unsigned int version)
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(MultiScattBrdfDriver)
     & FP_NVP_(brdf_interface);
  boost::serialization::split_member(ar, *this, version);
}

template<class Archive>
void LidortBrdfDriver::save(Archive & UNUSED(a), const unsigned int UNUSED(version)) const
{
  // Nothing more to do
}
template<class Archive>
void LidortBrdfDriver::load(Archive & UNUSED(ar), const unsigned int UNUSED(version))
{
  // Can not place in MultiScattBrdfDriver since brdf_interface_ is initalized here
  Spurr_Brdf_Sup_Inputs_Base& brdf_inputs = brdf_interface()->brdf_sup_inputs_base();
  brdf_params.reference( brdf_inputs.bs_brdf_parameters() );
  brdf_factors.reference( brdf_inputs.bs_brdf_factors() );
}

FP_IMPLEMENT(LidortBrdfDriver);

#endif

//-----------------------------------------------------------------------
/// Initialize LIDORT BRDF interface
//-----------------------------------------------------------------------

LidortBrdfDriver::LidortBrdfDriver(int nstream, int nmoment)
  : MultiScattBrdfDriver(nstream, nmoment)
{
  brdf_interface_.reset( new Brdf_Lin_Sup_Masters() );

  init();
}
