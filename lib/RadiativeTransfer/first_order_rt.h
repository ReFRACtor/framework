#ifndef FIRST_ORDER_RT_H
#define FIRST_ORDER_RT_H

#include "first_order_driver.h"
#include "spurr_rt.h"

// Include BRDF driver interface from LIDORT driver
#include "lidort_driver.h"

namespace FullPhysics {

/****************************************************************//**
  Uses the Spurr interfaces to construct a radiative transfer
  class connecting L2 FP and First Order
 *******************************************************************/
class FirstOrderRt : public SpurrRt {
public:
  FirstOrderRt(const boost::shared_ptr<RtAtmosphere>& Atm,
               const boost::shared_ptr<StokesCoefficient>& Stokes_coef,
               const blitz::Array<double, 1>& Sza, 
               const blitz::Array<double, 1>& Zen, 
               const blitz::Array<double, 1>& Azm,
               int Number_streams,
               int Number_moments,
               bool do_solar = true,
               bool do_thermal = false);

  /// Number of quadtature streams in the cosine half space
  virtual int number_stream() const { return rt_driver()->number_stream(); }
  
  /// Number of moments for scattering matrix.
  int number_moment() const { return rt_driver()->number_moment(); };
  
  /// Convenience routine to get brdf driver object
  const boost::shared_ptr<LidortBrdfDriver> brdf_driver() const
    { return boost::shared_ptr<LidortBrdfDriver>(boost::dynamic_pointer_cast<LidortBrdfDriver>(rt_driver()->brdf_driver())); }

  /// Convenience routine to get rt driver object
  const boost::shared_ptr<FirstOrderDriver> rt_driver() const 
  { return boost::shared_ptr<FirstOrderDriver>(boost::dynamic_pointer_cast<FirstOrderDriver>(rt_driver_)); }

  virtual void print(std::ostream& Os, bool Short_form = false) const;
private:
  FirstOrderRt() {}
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
}

FP_EXPORT_KEY(FirstOrderRt);
#endif
