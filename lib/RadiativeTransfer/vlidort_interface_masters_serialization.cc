#include "vlidort_interface_masters.h"
#include "fp_serialize_support.h"
#include "linear_algebra.h"

using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION

////////
// Serialization functions of VBrdf_Linsup_Masters
template<class Archive>
void VBrdf_Linsup_Masters::serialize(Archive& ar, const unsigned int version)
{
  FP_GENERIC_BASE(VBrdf_Linsup_Masters);
    ar 
       & FP_NVP_(vbrdf_sup_in)
       & FP_NVP_(vbrdf_linsup_in)
       & FP_NVP_(vbrdf_sup_inputstatus)
       & FP_NVP_(vbrdf_sup_out)
       & FP_NVP_(vbrdf_linsup_out)
       & FP_NVP_(vbrdf_sup_outputstatus);

  boost::serialization::split_member(ar, *this, version);
}

template<class Archive>
void VBrdf_Linsup_Masters::save(Archive & UNUSED(a), const unsigned int UNUSED(version)) const
{
  // Nothing more to do
}

template<class Archive>
void VBrdf_Linsup_Masters::load(Archive & UNUSED(ar), const unsigned int UNUSED(version))
{
  // The arrays come back from serialization in C order. We need them
  // in fortran order, so change them to that.
}

////////
// Serialization functions of VBrdf_Sup_Masters
template<class Archive>
void VBrdf_Sup_Masters::serialize(Archive& ar, const unsigned int version)
{
  FP_GENERIC_BASE(VBrdf_Sup_Masters);
    ar 
       & FP_NVP_(vbrdf_sup_in)
       & FP_NVP_(vbrdf_sup_inputstatus)
       & FP_NVP_(vbrdf_sup_out)
       & FP_NVP_(vbrdf_sup_outputstatus);

  boost::serialization::split_member(ar, *this, version);
}

template<class Archive>
void VBrdf_Sup_Masters::save(Archive & UNUSED(a), const unsigned int UNUSED(version)) const
{
  // Nothing more to do
}

template<class Archive>
void VBrdf_Sup_Masters::load(Archive & UNUSED(ar), const unsigned int UNUSED(version))
{
  // The arrays come back from serialization in C order. We need them
  // in fortran order, so change them to that.
}

////////
// Serialization functions of VLidort_Inputs
template<class Archive>
void VLidort_Inputs::serialize(Archive& ar, const unsigned int version)
{
  FP_GENERIC_BASE(VLidort_Inputs);
    ar 
       & FP_NVP_(vlidort_sup)
       & FP_NVP_(vlidort_fixin)
       & FP_NVP_(vlidort_modin)
       & FP_NVP_(vlidort_inputstatus);

  boost::serialization::split_member(ar, *this, version);
}

template<class Archive>
void VLidort_Inputs::save(Archive & UNUSED(a), const unsigned int UNUSED(version)) const
{
  // Nothing more to do
}

template<class Archive>
void VLidort_Inputs::load(Archive & UNUSED(ar), const unsigned int UNUSED(version))
{
  // The arrays come back from serialization in C order. We need them
  // in fortran order, so change them to that.
}

////////
// Serialization functions of VLidort_Masters
template<class Archive>
void VLidort_Masters::serialize(Archive& ar, const unsigned int version)
{
  FP_GENERIC_BASE(VLidort_Masters);
    ar 
       & FP_NVP_(vlidort_fixin)
       & FP_NVP_(vlidort_modin)
       & FP_NVP_(vlidort_sup)
       & FP_NVP_(vlidort_out);

  boost::serialization::split_member(ar, *this, version);
}

template<class Archive>
void VLidort_Masters::save(Archive & UNUSED(a), const unsigned int UNUSED(version)) const
{
  // Nothing more to do
}

template<class Archive>
void VLidort_Masters::load(Archive & UNUSED(ar), const unsigned int UNUSED(version))
{
  // The arrays come back from serialization in C order. We need them
  // in fortran order, so change them to that.
}

////////
// Serialization functions of VLidort_L_Inputs
template<class Archive>
void VLidort_L_Inputs::serialize(Archive& ar, const unsigned int version)
{
  FP_GENERIC_BASE(VLidort_L_Inputs);
    ar 
       & FP_NVP_(vlidort_linsup)
       & FP_NVP_(vlidort_fixin)
       & FP_NVP_(vlidort_modin)
       & FP_NVP_(vlidort_linfixin)
       & FP_NVP_(vlidort_linmodin)
       & FP_NVP_(vlidort_inputstatus);

  boost::serialization::split_member(ar, *this, version);
}

template<class Archive>
void VLidort_L_Inputs::save(Archive & UNUSED(a), const unsigned int UNUSED(version)) const
{
  // Nothing more to do
}

template<class Archive>
void VLidort_L_Inputs::load(Archive & UNUSED(ar), const unsigned int UNUSED(version))
{
  // The arrays come back from serialization in C order. We need them
  // in fortran order, so change them to that.
}

////////
// Serialization functions of VLidort_Lcs_Masters
template<class Archive>
void VLidort_Lcs_Masters::serialize(Archive& ar, const unsigned int version)
{
  FP_GENERIC_BASE(VLidort_Lcs_Masters);
    ar 
       & FP_NVP_(vlidort_fixin)
       & FP_NVP_(vlidort_modin)
       & FP_NVP_(vlidort_sup)
       & FP_NVP_(vlidort_out)
       & FP_NVP_(vlidort_linfixin)
       & FP_NVP_(vlidort_linmodin)
       & FP_NVP_(vlidort_linsup)
       & FP_NVP_(vlidort_linout);

  boost::serialization::split_member(ar, *this, version);
}

template<class Archive>
void VLidort_Lcs_Masters::save(Archive & UNUSED(a), const unsigned int UNUSED(version)) const
{
  // Nothing more to do
}

template<class Archive>
void VLidort_Lcs_Masters::load(Archive & UNUSED(ar), const unsigned int UNUSED(version))
{
  // The arrays come back from serialization in C order. We need them
  // in fortran order, so change them to that.
}

////////
// Serialization functions of VLidort_Lps_Masters
template<class Archive>
void VLidort_Lps_Masters::serialize(Archive& ar, const unsigned int version)
{
  FP_GENERIC_BASE(VLidort_Lps_Masters);
    ar 
       & FP_NVP_(vlidort_fixin)
       & FP_NVP_(vlidort_modin)
       & FP_NVP_(vlidort_sup)
       & FP_NVP_(vlidort_out)
       & FP_NVP_(vlidort_linfixin)
       & FP_NVP_(vlidort_linmodin)
       & FP_NVP_(vlidort_linsup)
       & FP_NVP_(vlidort_linout);

  boost::serialization::split_member(ar, *this, version);
}

template<class Archive>
void VLidort_Lps_Masters::save(Archive & UNUSED(a), const unsigned int UNUSED(version)) const
{
  // Nothing more to do
}

template<class Archive>
void VLidort_Lps_Masters::load(Archive & UNUSED(ar), const unsigned int UNUSED(version))
{
  // The arrays come back from serialization in C order. We need them
  // in fortran order, so change them to that.
}

////////
// Serialization functions of VLidort_Vbrdf_Sup_Accessories
template<class Archive>
void VLidort_Vbrdf_Sup_Accessories::serialize(Archive& ar, const unsigned int version)
{
  FP_GENERIC_BASE(VLidort_Vbrdf_Sup_Accessories);
    ar 
       & FP_NVP_(vbrdf_sup_out)
       & FP_NVP_(vlidort_fixin)
       & FP_NVP_(vlidort_modin)
       & FP_NVP_(vlidort_sup)
       & FP_NVP_(vbrdf_sup_in)
       & FP_NVP_(vlidort_vbrdfcheck_status);

  boost::serialization::split_member(ar, *this, version);
}

template<class Archive>
void VLidort_Vbrdf_Sup_Accessories::save(Archive & UNUSED(a), const unsigned int UNUSED(version)) const
{
  // Nothing more to do
}

template<class Archive>
void VLidort_Vbrdf_Sup_Accessories::load(Archive & UNUSED(ar), const unsigned int UNUSED(version))
{
  // The arrays come back from serialization in C order. We need them
  // in fortran order, so change them to that.
}

FP_IMPLEMENT(VBrdf_Linsup_Masters);
FP_IMPLEMENT(VBrdf_Sup_Masters);
FP_IMPLEMENT(VLidort_Inputs);
FP_IMPLEMENT(VLidort_Masters);
FP_IMPLEMENT(VLidort_L_Inputs);
FP_IMPLEMENT(VLidort_Lcs_Masters);
FP_IMPLEMENT(VLidort_Lps_Masters);
FP_IMPLEMENT(VLidort_Vbrdf_Sup_Accessories);

#endif