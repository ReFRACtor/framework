#ifndef AEROSOL_PROPERTY_IMP_BASE_H
#define AEROSOL_PROPERTY_IMP_BASE_H

#include "aerosol_property.h"
#include "sub_state_vector_array.h"
#include <boost/lexical_cast.hpp>

namespace FullPhysics {
/****************************************************************//**
  As a design principle, we have each base class with the absolutely
  minimum interface needed for use from the rest of the system. This
  allows us to support any future code that supports this minimum 
  interface.
  
  However, almost always you will want to derive from this class 
  instead. See PressureImpBase for a more complete discussion of this.
*******************************************************************/
class AerosolPropertyImpBase:
    virtual public SubStateVectorArray<AerosolProperty>
{
public:
  virtual ~AerosolPropertyImpBase() {}
  virtual boost::shared_ptr<AerosolProperty> clone() const = 0;
  virtual ArrayAd<double,1> extinction_coefficient_each_layer(double wn) 
    const = 0;
  virtual ArrayAd<double, 1> scattering_coefficient_each_layer(double wn)
    const = 0;
  virtual ArrayAd<double, 3> 
  phase_function_moment_each_layer(double wn, int nmom = -1, 
				   int nscatt = -1) const = 0;

  virtual std::string sub_state_identifier() const { return "aerosol/property"; }

  virtual std::string state_vector_name_i(int i) const 
  {
      return "Aerosol Property Coeff " + boost::lexical_cast<std::string>(i + 1);
  }
  virtual void print(std::ostream& Os) const
  { Os << "AerosolPropertyImpBase"; }

//-----------------------------------------------------------------------
/// Returns the value of the coefficients used to generate the aerosol
/// property
//-----------------------------------------------------------------------

  blitz::Array<double, 1> aerosol_parameter() const
  {
    return coefficient().value();
  }

//-----------------------------------------------------------------------
/// Returns the uncertainty of the aerosol type coefficients
//-----------------------------------------------------------------------

  blitz::Array<double, 1> aerosol_parameter_uncertainty() const
  {
    blitz::Array<double, 1> uncert(coefficient().rows());
    for(int i = 0; i < sv_cov_sub.rows(); i++)
      uncert(i) = (sv_cov_sub(i,i) > 0 ? sqrt(sv_cov_sub(i, i)) : 0.0);
    return uncert;
  }

protected:
  // Don't think we need a cache, but if we end up needing it can add
  // this like we have in AerosolExtinctionImpBase.
  // mutable bool cache_stale;

//-----------------------------------------------------------------------
/// Initialize object.
//-----------------------------------------------------------------------

  void init(const blitz::Array<double, 1>& Coeff)
  { 
    SubStateVectorArray<AerosolProperty>::init(Coeff);
  }

//-----------------------------------------------------------------------
/// Default constructor, derived class should call init if they use this
/// constructor.
//-----------------------------------------------------------------------

  AerosolPropertyImpBase() { }

//-----------------------------------------------------------------------
/// Constructor that sets the coefficient() values.
//-----------------------------------------------------------------------
  AerosolPropertyImpBase(const blitz::Array<double, 1>& Coeff)
  {
    SubStateVectorArray<AerosolProperty>::init(Coeff);
  }
private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
typedef SubStateVectorArray<AerosolProperty> SubStateVectorArrayAerosolProperty;
}

FP_EXPORT_KEY(AerosolPropertyImpBase);
FP_EXPORT_KEY(SubStateVectorArrayAerosolProperty);

#endif
