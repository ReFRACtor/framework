#ifndef ABSORBER_VMR_LEVEL_H
#define ABSORBER_VMR_LEVEL_H

#include "absorber_vmr_imp_base.h"

namespace FullPhysics {
/****************************************************************//**
  This class maps the state vector to the absorber VMR on each
  level.

  This particular implementation uses the state vector values as
  the VMR for each pressure level.
*******************************************************************/
class AbsorberVmrLevel : virtual public AbsorberVmrImpBase {
public:
  AbsorberVmrLevel(const boost::shared_ptr<Pressure>& Mapped_Press,
		   const blitz::Array<double, 1>& Vmr,
		   const std::string& Gas_name,
		   boost::shared_ptr<StateMapping> in_map = boost::make_shared<StateMappingLinear>(),
		   const boost::shared_ptr<Pressure>& Coeff_Press = NULL);

  virtual ~AbsorberVmrLevel() {}
  virtual void print(std::ostream& Os) const;
  virtual std::string sub_state_identifier() const
  {
    return "absorber_levels/" + mapping->name() + "/" + gas_name();
  }

  virtual std::string state_vector_name_i(int i) const;

  virtual boost::shared_ptr<AbsorberVmr> clone() const;

  //-----------------------------------------------------------------------
  /// Pressure levels that the mapped vmr is on
  //-----------------------------------------------------------------------
  blitz::Array<double, 1> pressure_profile() const;
  
  //-----------------------------------------------------------------------
  /// VMR on the mapped pressure grid.
  //-----------------------------------------------------------------------
  blitz::Array<double, 1> vmr_profile() const;
  
  //-----------------------------------------------------------------------
  /// Pressure levels that the retrieval vmr cofficients are on
  //-----------------------------------------------------------------------
  blitz::Array<double, 1> coeff_pressure_profile() const;
  
protected:
  virtual void calc_vmr() const;
  AbsorberVmrLevel() {}
private:
  boost::shared_ptr<Pressure> coeff_pressure;
  
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
}

FP_EXPORT_KEY(AbsorberVmrLevel);
#endif
