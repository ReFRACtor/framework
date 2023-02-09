#include "vlidort_brdf_driver.h"

#include "fp_serialize_support.h"
#include "fp_exception.h"
#include "linear_algebra.h"
#include "ostream_pad.h"
#include "fe_disable_exception.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef FP_HAVE_BOOST_SERIALIZATION

template<class Archive>
void VLidortBrdfDriver::serialize(Archive & ar,
                                  const unsigned int version)
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(MultiScattBrdfDriver)
    & FP_NVP_(brdf_interface);
  boost::serialization::split_member(ar, *this, version);
}

template<class Archive>
void VLidortBrdfDriver::save(Archive & UNUSED(a), const unsigned int UNUSED(version)) const
{
  // Nothing more to do
}
template<class Archive>
void VLidortBrdfDriver::load(Archive & UNUSED(ar), const unsigned int UNUSED(version))
{
  // Can not place in MultiScattBrdfDriver since brdf_interface_ is initalized here
  VBrdf_Sup_Inputs& brdf_inputs = brdf_interface_->vbrdf_sup_in();
  brdf_params.reference( brdf_inputs.bs_brdf_parameters() );
  brdf_factors.reference( brdf_inputs.bs_brdf_factors() );
}

FP_IMPLEMENT(VLidortBrdfDriver);

#endif

//-----------------------------------------------------------------------
/// Initialize VLIDORT BRDF interface
//-----------------------------------------------------------------------

VLidortBrdfDriver::VLidortBrdfDriver(int nstream, int nmoment)
  : MultiScattBrdfDriver(nstream, nmoment)
{
  brdf_interface_.reset( new VBrdf_Linsup_Masters() );

  init();
}
