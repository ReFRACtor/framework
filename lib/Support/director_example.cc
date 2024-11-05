#include "director_example.h"
#include "fp_serialize_support.h"

using namespace FullPhysics;

template<class Archive>
void DirectorExample::serialize(Archive & ar, const unsigned int )
{
  FP_GENERIC_BASE(DirectorExample);
  ar & FP_NVP_(value_to_add) & FP_NVP_(name);
  // Note, this is just here so we can have a weak_ptr to a
  // DirectorExample handled correctly. If you are looking at this as
  // an example, you probably want to do something a little higher
  // order like the Observer/Observable pattern in ReFRACtor rather
  // than directly calling add_ptr_serialized_reference. But for this
  // skeleton of swig, we don't want to put all that code in place so
  // we just do a low level example.
  add_ptr_serialized_reference(this);
}

template<class Archive>
void DirectorExampleUser::serialize(Archive & ar, const unsigned int )
{
  FP_GENERIC_BASE(DirectorExampleUser);
  ar & FP_NVP_(value_vec) & FP_NVP_(value_vec2)
    & FP_NVP_(value) & FP_NVP_(t) & FP_NVP_(t2);
}

template<class Archive>
void DirectorExampleWeakPtr::serialize(Archive & ar, const unsigned int )
{
  FP_GENERIC_BASE(DirectorExampleWeakPtr);
  ar & FP_NVP_(a);
}

FP_IMPLEMENT(DirectorExample);
FP_IMPLEMENT(DirectorExampleUser);
FP_IMPLEMENT(DirectorExampleWeakPtr);

void DirectorExample::print(std::ostream& Os) const
{
  Os << "DirectorExample:\n"
     << "  name:         " << name_
     << "  value_to_add: " << value_to_add_;
}
