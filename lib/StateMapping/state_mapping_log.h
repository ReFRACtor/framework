#ifndef STATE_MAPPING_LOG_H
#define STATE_MAPPING_LOG_H

#include <blitz/array.h>
#include <boost/shared_ptr.hpp>

#include "array_ad.h"
#include "state_mapping.h"


namespace FullPhysics {
/****************************************************************//**
  This class implements log encoding of coeffs for the retrieval view
  while using the real (linear) values for the forward model view.

  For additional information see docs for StateMapping class.
*******************************************************************/
class StateMappingLog : public StateMapping  {
public:
  //-----------------------------------------------------------------------
  /// Default Constructor.
  //-----------------------------------------------------------------------

  StateMappingLog() {};
  virtual ~StateMappingLog() {}
  
  virtual ArrayAd<double, 1> mapped_state(const ArrayAd<double, 1>& retrieval_values)
    const {
    blitz::Array<AutoDerivative<double>, 1> res(retrieval_values.rows());
    for(int i = 0; i < res.rows(); ++i)
      res(i) = std::exp(retrieval_values(i));
    return ArrayAd<double, 1>(res);
  }
  
  virtual ArrayAd<double, 1> retrieval_state
  (const ArrayAd<double, 1>& initial_values) const {
    blitz::Array<AutoDerivative<double>, 1> res(initial_values.rows());
    for(int i = 0; i < res.rows(); ++i)
      res(i) = std::log(initial_values(i));
    return ArrayAd<double, 1>(res);
  }

  virtual std::string name() const { return "log"; }

  virtual boost::shared_ptr<StateMapping> clone() const
  { return boost::shared_ptr<StateMapping>(new StateMappingLog()); }

private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
}

FP_EXPORT_KEY(StateMappingLog);

#endif
