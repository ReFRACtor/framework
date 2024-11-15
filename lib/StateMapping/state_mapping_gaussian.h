#ifndef STATE_MAPPING_GAUSSIAN_H
#define STATE_MAPPING_GAUSSIAN_H

#include <blitz/array.h>

#include "array_ad.h"
#include "state_mapping.h"
#include "pressure.h"

namespace FullPhysics {
/****************************************************************//**
  This class implements Gaussian parameterization of coeffs for the
  retrieval view  while using the calculated values for the forward
  model view.

  Note that due to the scaling by total optical depth performed in
  this class, it is only suitable for Absorbers.

  For additional information see docs for StateMapping class.
*******************************************************************/
class StateMappingGaussian : public StateMapping  {
public:
  StateMappingGaussian(const boost::shared_ptr<Pressure>& in_press,
                       bool Linear_Total,
                       double Min_Desired = 1e-9);
  virtual ~StateMappingGaussian() {};

  virtual boost::shared_ptr<StateMapping> clone() const
  {
    return boost::shared_ptr<StateMapping>(new StateMappingGaussian(press, linear_total));
  }

  //-----------------------------------------------------------------------
  /// Whether this mapping uses a linear total parameter (alternative is log)
  //-----------------------------------------------------------------------

  virtual bool is_linear_total() const { return linear_total; }
  double min_desired() const {return min_desired_;}
  const boost::shared_ptr<Pressure>& pressure() const { return press; }
  
  virtual AutoDerivative<double> total_optical_depth
  (ArrayAd<double, 1> component) const;

  virtual ArrayAd<double, 1> mapped_state
  (const ArrayAd<double, 1>& retrieval_values) const;
  
  virtual ArrayAd<double, 1> retrieval_state
  (const ArrayAd<double, 1>& initial_values) const
  { return initial_values;}

  virtual std::string name() const { return map_name; }
private:
  std::string map_name;
  double min_desired_;
  bool linear_total;

  //-----------------------------------------------------------------------
  /// Pressure portion of the state
  /// Note that levels that define the layers used in the Radiative Transfer
  /// calculation may vary as we do a retrieval.
  //-----------------------------------------------------------------------
  boost::shared_ptr<Pressure> press;
  StateMappingGaussian() {}
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
}

FP_EXPORT_KEY(StateMappingGaussian);
#endif
