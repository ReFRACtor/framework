#ifndef GROUND_LAMBERTIAN_PIECEWISE_H
#define GROUND_LAMBERTIAN_PIECEWISE_H

#include "ground_piecewise.h"

namespace FullPhysics {

/****************************************************************//**
  This class implements an lambertian implemented as a piecewise
  linear interpolation. This would be a single value that is different
  for each spectral point. 
*******************************************************************/
class GroundLambertianPiecewise: virtual public GroundPiecewise {

public:
  GroundLambertianPiecewise(const ArrayWithUnit<double, 1>& spectral_points,
                            const blitz::Array<double, 1>& point_values);

  virtual boost::shared_ptr<Ground> clone() const;

  virtual std::string sub_state_identifier() const;
  virtual std::string state_vector_name_i(int i) const;
  virtual void print(std::ostream& Os) const;
  virtual std::string desc() const;
private:
  GroundLambertianPiecewise() {}
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
}

FP_EXPORT_KEY(GroundLambertianPiecewise);
#endif
