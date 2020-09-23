#ifndef STATE_MAPPING_LINEAR_H
#define STATE_MAPPING_LINEAR_H

#include "state_mapping.h"
#include "array_ad.h"

namespace FullPhysics {
/****************************************************************//**
  This class implements linear 1-to-1 mapping and just provides the
  SubStateVectorArray's coefficients as-is.

  For additional information see docs for StateMapping class.
*******************************************************************/
class StateMappingLinear : public StateMapping {
public:
  virtual ~StateMappingLinear() {}

  //-----------------------------------------------------------------------
  /// Default Constructor.
  //-----------------------------------------------------------------------

  StateMappingLinear() {}

  //-----------------------------------------------------------------------
  /// Calculation of forward model view of coeffs with mapping applied
  //-----------------------------------------------------------------------

  virtual ArrayAd<double, 1> mapped_state
  (const ArrayAd<double, 1>& retrieval_values) const
  { return retrieval_values; }

  //-----------------------------------------------------------------------
  /// Calculation of initial retrieval view  of coeffs with mapping applied
  //-----------------------------------------------------------------------

  virtual ArrayAd<double, 1> retrieval_state
  (const ArrayAd<double, 1>& initial_values) const
  { return initial_values; }

  //-----------------------------------------------------------------------
  /// Index into initial values for each retrieval state entry
  //-----------------------------------------------------------------------
  virtual int initial_values_index(const int retrieval_state_index) const
  { return retrieval_state_index; }

  //-----------------------------------------------------------------------
  /// Assigned mapping name
  //-----------------------------------------------------------------------
  
  virtual std::string name() const { return "linear"; }

  virtual boost::shared_ptr<StateMapping> clone() const
  {
    return boost::shared_ptr<StateMapping>(new StateMappingLinear());
  }
private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
}

FP_EXPORT_KEY(StateMappingLinear);
#endif
