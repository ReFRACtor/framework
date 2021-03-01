#ifndef SAMPLE_GRID_IMP_BASE_H
#define SAMPLE_GRID_BASE_H
#include "sample_grid.h"
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
class SampleGridImpBase: virtual public SubStateVectorArray<SampleGrid> {
public:
  virtual ~SampleGridImpBase() = default;
  virtual boost::shared_ptr<SampleGrid> clone() const = 0;
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
  virtual std::string desc() const { return "SampleGridImpBase"; }

  SampleGridImpBase(const SampleGridImpBase& V) = default;

protected:

//-----------------------------------------------------------------------
/// Default constructor, derived class should call init if they use this
/// constructor.
//-----------------------------------------------------------------------

  SampleGridImpBase() { }

//-----------------------------------------------------------------------
/// Constructor that sets the coefficient() values.
/// See SubStateVectorArray for a discussion of Mark_according_to_press and
/// Pdep_start.
//-----------------------------------------------------------------------
  SampleGridImpBase(const blitz::Array<double, 1>& Coeff, 
                    boost::shared_ptr<StateMapping> in_map = boost::make_shared<StateMappingLinear>())
  {
    init(Coeff, in_map);
  }

private:

  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);

};
typedef SubStateVectorArray<SampleGrid> SubStateVectorArraySampleGrid;
}

FP_EXPORT_KEY(SampleGridImpBase);
FP_EXPORT_KEY(SubStateVectorArraySampleGrid);
#endif
