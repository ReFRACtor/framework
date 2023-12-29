#ifndef GROUND_PIECEWISE_H
#define GROUND_PIECEWISE_H

#include "ground_imp_base.h"
#include "array_with_unit.h"
#include "linear_interpolate.h"
#include "state_mapping_linear.h"

namespace FullPhysics {

/****************************************************************//**
  This class implements a ground type implemented as a piecewise
  linear interpolation. This would be a single value that is different
  for each spectral point. Subclasses should add the state vector
  identifcation information for the type of ground information it
  represents. This class is absctract.
*******************************************************************/
class GroundPiecewise: virtual public GroundImpBase {
public:
  GroundPiecewise(const ArrayWithUnit<double, 1>& spectral_points,
                  const blitz::Array<double, 1>& point_values,
                  const boost::shared_ptr<StateMapping>& mapping = boost::make_shared<StateMappingLinear>());

  virtual const ArrayWithUnit<double, 1>& spectral_points() const;

  virtual ArrayAd<double, 1> surface_parameter(const double wn, const int spec_index) const;

  virtual const AutoDerivative<double> value_at_point(const DoubleWithUnit wave_point) const;

  virtual void update_sub_state_hook();

  virtual boost::shared_ptr<Ground> clone() const = 0;

  virtual std::string sub_state_identifier() const = 0;
  virtual std::string state_vector_name_i(int i) const = 0;

  virtual void print(std::ostream& Os) const {Os << "GroundPiecewise";}
  GroundPiecewise(const GroundPiecewise& V) = default;
protected:
  ArrayWithUnit<double, 1> spectral_points_;
  GroundPiecewise() {};
private:
  // Interpolation object is update for each state vector update
  boost::shared_ptr<LinearInterpolate<AutoDerivative<double>, AutoDerivative<double> > > ground_interp;
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
  template<class Archive>
  void save(Archive & ar, const unsigned int version) const;
  template<class Archive>
  void load(Archive & ar, const unsigned int version);
};
}

FP_EXPORT_KEY(GroundPiecewise);
#endif
