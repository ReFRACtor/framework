%include "fp_common.i"
%{
#include "aerosol_extinction_imp_base.h"
#include "sub_state_vector_array.h"
%}

%base_import(aerosol_extinction)
%base_import(sub_state_vector_array)
%base_import(state_mapping)
%base_import(state_mapping_linear)

%fp_shared_ptr(FullPhysics::AerosolExtinctionImpBase);
%fp_shared_ptr(FullPhysics::SubStateVectorArray<FullPhysics::AerosolExtinction>);
namespace FullPhysics {
%template(SubStateVectorArrayAerosolExtinction) 
     FullPhysics::SubStateVectorArray<AerosolExtinction>;

// Allow these classes to be derived from in Python.
%feature("director") AerosolExtinctionImpBase;

// Note, at least for SWIG 2.0.4 a class that is derived from in python 
// needs to declare every virtual function that can be called on it, even 
// if all that happens is the base class to a director is called. This is 
// because this class is used to create the SwigDirector class, and this 
// class needs each of the member functions to direct things properly. It 
// is *not* necessary to add these function to the underlying
// C++, only that you declare them here.
//
// If you miss something, then you will get something like a recursion
// error in python when a virtual function is used that isn't explicitly
// listed here.
//
// This seems like a bug in 2.0.4, if SWIG needs all the member functions
// it should know to make them itself. So perhaps a future version of SWIG
// won't have this same constraint. But for now, this is required.

class AerosolExtinctionImpBase: public SubStateVectorArray<AerosolExtinction> {
public:
  // From AerosolExtinctionImpBase
  virtual ~AerosolExtinctionImpBase();
  virtual std::string desc() const;
  virtual boost::shared_ptr<AerosolExtinction> clone() const = 0;
  virtual AutoDerivative<double> extinction_for_layer(int i) const;
  %python_attribute(aerosol_name, virtual std::string);
  virtual ArrayAd<double, 1> aerosol_extinction(Pressure::PressureGridType Gtype = Pressure::INCREASING_PRESSURE) const;
  %python_attribute_abstract(model_short_name, std::string);
  %sub_state_virtual_func(AerosolExtinction);
  %python_attribute(aerosol_parameter, blitz::Array<double, 1>);
  %python_attribute(aerosol_parameter_uncertainty, blitz::Array<double, 1>);
  %pickle_serialization();
protected:
  mutable bool cache_stale;
  mutable ArrayAd<double, 1> aext;
  virtual void calc_aerosol_extinction() const = 0;
  %python_attribute(total_aod, AutoDerivative<double>);
  AerosolExtinctionImpBase(const std::string& Aerosol_name,
                           const blitz::Array<double, 1>& Coeff,
                           const boost::shared_ptr<Pressure>& Press,
                           boost::shared_ptr<StateMapping> in_map = boost::make_shared<StateMappingLinear>());
};
}

%template(vector_aerosol_extinction_imp_base) std::vector<boost::shared_ptr<FullPhysics::AerosolExtinctionImpBase> >;

// Extra code for handling boost serialization/python pickle of
// director classes
%{
// Needed by code below, can't easily figure these names out
// automatically so just include here
#include "aerosol_extinction_imp_base_wrap.h"
%}
%fp_director_serialization(AerosolExtinctionImpBase)

// List of things "import *" will include
%python_export("AerosolExtinctionImpBase");
