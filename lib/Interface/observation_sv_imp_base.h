#ifndef OBSERVATION_SV_IMP_BASE_H
#define OBSERVATION_SV_IMP_BASE_H
#include "observation_sv.h"
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
class ObservationSvImpBase: virtual public SubStateVectorArray<ObservationSv> {
public:
  virtual ~ObservationSvImpBase() {}
  virtual boost::shared_ptr<ObservationSv> clone() const = 0;
  virtual void update_sub_state_hook() 
  { // Nothing to do
  }
  
//-----------------------------------------------------------------------
/// Print to stream. The default calls the function "desc" that returns
/// a string. This gives cleaner interface for deriving from this class
/// in python, but most C++ classes will want to override this function
/// rather than using desc.
//-----------------------------------------------------------------------
  virtual void print(std::ostream& Os) const { Os << desc(); }

//-----------------------------------------------------------------------
/// Description of object, to be printed to stream. This gives a cleaner
/// interface for deriving from python.
//-----------------------------------------------------------------------
  virtual std::string desc() const { return "ObservationSvImpBase"; }

  virtual blitz::Array<bool, 1> state_used() const
  {
    blitz::Array<bool, 1> res(sv_full.rows());
    res = false;
    if(state_vector_start_index() < res.rows())
      res(blitz::Range(state_vector_start_index(), 
           state_vector_start_index() + sub_vector_size() - 1)) = true;
    return res;
  }
  ObservationSvImpBase(const ObservationSvImpBase& V) = default;
protected:
//-----------------------------------------------------------------------
/// Default constructor, derived class should call init if they use this
/// constructor.
//-----------------------------------------------------------------------

  ObservationSvImpBase() { }

//-----------------------------------------------------------------------
/// Constructor that sets the coefficient() values.
/// See SubStateVectorArray for a discussion of Mark_according_to_press and
/// Pdep_start.
//-----------------------------------------------------------------------
  ObservationSvImpBase(const blitz::Array<double, 1>& Coeff, 
                boost::shared_ptr<StateMapping> in_map =
                  boost::make_shared<StateMappingLinear>())
  {
    init(Coeff, in_map);
  }
private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
typedef SubStateVectorArray<ObservationSv> SubStateVectorArrayObservationSv;
}

FP_EXPORT_KEY(ObservationSvImpBase);
FP_EXPORT_KEY(SubStateVectorArrayObservationSv);
#endif
