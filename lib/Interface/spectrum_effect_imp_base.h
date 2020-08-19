#ifndef SPECTRUM_EFFECT_IMP_BASE_H
#define SPECTRUM_EFFECT_IMP_BASE_H
#include "spectrum_effect.h"
#include "sub_state_vector_array.h"

namespace FullPhysics {
/****************************************************************//**
  As a design principle, we have each base class with the absolutely
  minimum interface needed for use from the rest of the system. This
  allows us to support any future code that supports this minimum 
  interface.
  
  However, almost always you will want to derive from this class 
  instead. See PressureImpBase for a more complete discussion of this.
*******************************************************************/
class SpectrumEffectImpBase:
    virtual public SubStateVectorArray<SpectrumEffect> {
public:
  virtual ~SpectrumEffectImpBase() {}
  virtual boost::shared_ptr<SpectrumEffect> clone() const = 0;

//-----------------------------------------------------------------------
/// Print to stream. The default calls the function "desc" that returns
/// a string. This gives cleaner interface for deriving from this class
/// in python, but most C++ classes will want to override this function
/// rather than using desc.
//-----------------------------------------------------------------------
  virtual void print(std::ostream& Os, bool UNUSED(Short_form) = false) const { Os << desc(); }

//-----------------------------------------------------------------------
/// Description of object, to be printed to stream. This gives a cleaner
/// interface for deriving from python.
//-----------------------------------------------------------------------
  virtual std::string desc() const { return "SpectrumEffectImpBase"; }

protected:

//-----------------------------------------------------------------------
/// Default constructor, derived class should call init if they use this
/// constructor.
//-----------------------------------------------------------------------
  SpectrumEffectImpBase() {}

//-----------------------------------------------------------------------
/// Constructor that sets the coefficient() values.
//-----------------------------------------------------------------------
  SpectrumEffectImpBase(const blitz::Array<double, 1>& Coeff, boost::shared_ptr<StateMapping> Mapping)
  {
    init(Coeff, Mapping);
  }

//-----------------------------------------------------------------------
/// Constructor that sets the coefficient values when
/// they are scalar values
//-----------------------------------------------------------------------
  SpectrumEffectImpBase(double Coeff)
  { init(Coeff); }
private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
typedef SubStateVectorArray<SpectrumEffect> SubStateVectorArraySpectrumEffect;
}

FP_EXPORT_KEY(SpectrumEffectImpBase);
FP_EXPORT_KEY(SubStateVectorArraySpectrumEffect);
#endif
