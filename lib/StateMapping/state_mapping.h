#ifndef STATE_MAPPING_H
#define STATE_MAPPING_H
#include "array_ad.h"
#include "pressure.h"
#include <string>

namespace FullPhysics {
/****************************************************************//**
  This class manages mapping between the State Vector and forward
  model.

  Individual components of the State Vector may have different
  representations (e.g. to add a scaling factor or log encode the
  values) when being retrieved/perturbed vs. the representations
  used in forward model calculations.

  Derived classes capture those different representations.
*******************************************************************/
class StateMapping : public Printable<StateMapping> {
public:
  StateMapping() {}
  virtual ~StateMapping() {};

  //-----------------------------------------------------------------------
  /// Calculation of forward model view of coeffs with mapping applied
  //-----------------------------------------------------------------------
  virtual ArrayAd<double, 1> mapped_state(const ArrayAd<double, 1>& retrieval_values) const = 0;

  //-----------------------------------------------------------------------
  /// Calculation of initial retrieval view of coeffs with mapping applied
  //-----------------------------------------------------------------------

  virtual ArrayAd<double, 1> retrieval_state(const ArrayAd<double, 1>& initial_values) const = 0;

  //-----------------------------------------------------------------------
  /// Index into initial values for each retrieval state entry if applicable
  //-----------------------------------------------------------------------
  virtual int initial_values_index(const int retrieval_state_index) const = 0;

  //-----------------------------------------------------------------------
  /// Assigned mapping name
  //-----------------------------------------------------------------------

  virtual std::string name() const = 0;

  virtual boost::shared_ptr<StateMapping> clone() const = 0;

  void print(std::ostream& Os) const
  { Os << "StateMapping\n"
       << "  name: " << name() << "\n";
  }
private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
}

FP_EXPORT_KEY(StateMapping);
#endif
