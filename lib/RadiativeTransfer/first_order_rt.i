%include "fp_common.i"
%{
#include "first_order_rt.h"
#include "sub_state_vector_array.h"
#include "pressure.h"
%}
%base_import(observer)
%base_import(rt_atmosphere)
%base_import(spurr_rt)
%import "lidort_driver.i"
%import "first_order_driver.i"
%fp_shared_ptr(FullPhysics::FirstOrderRt);

namespace FullPhysics {

// Force to be not abstract
%feature("notabstract") FirstOrderRt;

class FirstOrderRt : public SpurrRt  {
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


  %python_attribute(number_stream, virtual int)
  %python_attribute(number_moment, int)
  
  %python_attribute(brdf_driver, boost::shared_ptr<LidortBrdfDriver>)
  %python_attribute(rt_driver, boost::shared_ptr<FirstOrderDriver>)
  
  virtual void print(std::ostream& Os, bool Short_form = false) const;
  %pickle_serialization();
};
}
