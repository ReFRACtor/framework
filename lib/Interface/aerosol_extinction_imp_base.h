#ifndef AEROSOL_EXTINCTION_IMP_BASE_H
#define AEROSOL_EXTINCTION_IMP_BASE_H
#include "aerosol_extinction.h"
#include "sub_state_vector_array.h"
#include "state_mapping_linear.h"

namespace FullPhysics {
/****************************************************************//**
  As a design principle, we have each base class with the absolutely
  minimum interface needed for use from the rest of the system. This
  allows us to support any future code that supports this minimum 
  interface.
  
  However, almost always you will want to derive from this class 
  instead. See PressureImpBase for a more complete discussion of this.
*******************************************************************/
class AerosolExtinctionImpBase:
    virtual public SubStateVectorArray<AerosolExtinction> {
public:
  virtual ~AerosolExtinctionImpBase() {}
  virtual std::string aerosol_name() const {return aerosol_name_;}
  virtual ArrayAd<double, 1> aerosol_extinction() const
  { fill_cache(); return aext; }
  virtual AutoDerivative<double> extinction_for_layer(int i) const
  { fill_cache(); range_check(i, 0, aext.rows() - 1); 
    return (aext(i) + aext(i + 1)) / 2; }
  virtual boost::shared_ptr<AerosolExtinction> clone() const = 0;
  virtual void update_sub_state_hook() 
  { cache_stale = true; }

//-----------------------------------------------------------------------
/// A short name representing the type of extinction model being 
/// implemented
//-----------------------------------------------------------------------
  virtual std::string model_short_name() const = 0;
  
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
  virtual std::string desc() const { return "AerosolExtinctionImpBase"; }

//-----------------------------------------------------------------------
/// Returns the value of the coefficients used to generate the aerosol
/// extinction
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
//-----------------------------------------------------------------------
/// If this is true, the recalculate the vmr the next time we
/// need it.
//-----------------------------------------------------------------------
  mutable bool cache_stale;

//-----------------------------------------------------------------------
/// The cached aerosol extinction for each level. This should be
/// filled in by derived classes when calc_aerosol_extinction() is called.
//-----------------------------------------------------------------------
  mutable ArrayAd<double, 1> aext;

//-----------------------------------------------------------------------
/// Derived classes should provide a function to fill in vmr when this is 
/// called.
//-----------------------------------------------------------------------
  virtual void calc_aerosol_extinction() const = 0;

//-----------------------------------------------------------------------
/// Total aerosol optical depth of the extinction values in aext.
//-----------------------------------------------------------------------
  virtual AutoDerivative<double> total_aod() const
  {
    ArrayAd<double, 1> pressure_grid = press->pressure_grid().value;
    AutoDerivative<double> tot_aod = 0.0;
    for(int layer_idx = 0; layer_idx < press->number_layer(); layer_idx++) {
      AutoDerivative<double> delta_press = (pressure_grid(layer_idx + 1) - pressure_grid(layer_idx)) / 2.0;
      tot_aod += (delta_press * (aext(layer_idx) + aext(layer_idx + 1) ));
    }
    return tot_aod;
  }

//-----------------------------------------------------------------------
/// Initialize object.
//-----------------------------------------------------------------------

  void init(const std::string Aerosol_name,
            const blitz::Array<double, 1>& Coeff, 
            const boost::shared_ptr<Pressure>& Press,
            boost::shared_ptr<StateMapping> in_map = boost::make_shared<StateMappingLinear>())
  { 
    SubStateVectorArray<AerosolExtinction>::init(Coeff, in_map);
    aerosol_name_ = Aerosol_name;
    press = Press;
  }

//-----------------------------------------------------------------------
/// Default constructor, derived class should call init if they use this
/// constructor.
//-----------------------------------------------------------------------

  AerosolExtinctionImpBase() : cache_stale(true) { }

//-----------------------------------------------------------------------
/// Constructor that sets the coefficient() values.
//-----------------------------------------------------------------------
  AerosolExtinctionImpBase
  (const std::string& Aerosol_name,
   const blitz::Array<double, 1>& Coeff, 
   const boost::shared_ptr<Pressure>& Press,
   boost::shared_ptr<StateMapping> in_map = boost::make_shared<StateMappingLinear>())
    : cache_stale(true)
  {
    init(Aerosol_name, Coeff, Press, in_map);
  }
protected:
  boost::shared_ptr<Pressure> press;
private:
  void fill_cache() const
  {
    if(cache_stale) {
      aext.resize_number_variable(coeff.number_variable());
      aext.jacobian() = 0;
      calc_aerosol_extinction();
    }
    cache_stale = false;
  }
  std::string aerosol_name_;
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
typedef SubStateVectorArray<AerosolExtinction> SubStateVectorArrayAerosolExtinction;
}

FP_EXPORT_KEY(AerosolExtinctionImpBase);
FP_EXPORT_KEY(SubStateVectorArrayAerosolExtinction);
#endif
