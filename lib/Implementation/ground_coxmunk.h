#ifndef GROUND_COXMUNK_H
#define GROUND_COXMUNK_H

#include "ground_imp_base.h"
#include "auto_derivative.h"
#include "state_mapping_linear.h"

namespace FullPhysics {
/****************************************************************//**
  This class implements a Coxmunk ground type. 
*******************************************************************/
class GroundCoxmunk: virtual public GroundImpBase {
public:
  GroundCoxmunk(const double Windspeed,
                const blitz::Array<double, 1>& Refr_index,
                boost::shared_ptr<StateMapping> Mapping = boost::make_shared<StateMappingLinear>());

  virtual SpurrBrdfType spurr_brdf_type() const
  { return SpurrBrdfType::COXMUNK; }
  virtual ArrayAd<double, 1> surface_parameter(const double wn, const int spec_index) const;

  virtual const AutoDerivative<double> windspeed() const;

  virtual double refractive_index(const int Spec_idx) const;
  
  virtual boost::shared_ptr<Ground> clone() const;

  virtual std::string sub_state_identifier() const { return "ground/coxmunk"; }

  virtual std::string state_vector_name_i(int i) const;

  virtual void print(std::ostream& Os) const;
  
  virtual void update_sub_state_hook();

private:

  blitz::Array<double, 1> refractive_index_;
  GroundCoxmunk() {}
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
}

FP_EXPORT_KEY(GroundCoxmunk);
#endif
