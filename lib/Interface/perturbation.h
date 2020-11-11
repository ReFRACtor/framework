#ifndef PERTURBATION_H
#define PERTURBATION_H
#include "printable.h"
#include <blitz/array.h>

namespace FullPhysics {
/****************************************************************//**
  This gets the perturbation to use with a finite difference Jacobian.
*******************************************************************/
class Perturbation : public Printable<Perturbation> {
public:
  virtual ~Perturbation() {}

//-----------------------------------------------------------------------
/// Return the perturbation vector to use.
//-----------------------------------------------------------------------

  virtual blitz::Array<double, 1> perturbation() const = 0;

//-----------------------------------------------------------------------
/// Print description of object.
//-----------------------------------------------------------------------

  virtual void print(std::ostream& Os) const { Os << "Perturbation";}
private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
}

FP_EXPORT_KEY(Perturbation)
#endif
