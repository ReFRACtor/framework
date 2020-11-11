#ifndef SURFACE_TEMPERATURE_DIRECT_H
#define SURFACE_TEMPERATURE_DIRECT_H

#include "surface_temperature.h"
#include "sub_state_vector_array.h"
#include "double_with_unit.h"
#include "auto_derivative_with_unit.h"

namespace FullPhysics {
/****************************************************************//**
 Implements a direct representation of the surface temperature
 in the atmospheric state from the supplied value.
*******************************************************************/
class SurfaceTemperatureDirect :
    virtual public SubStateVectorArray<SurfaceTemperature> {
public:
  SurfaceTemperatureDirect(const ArrayWithUnit<double, 1>& surf_temp);
  virtual ~SurfaceTemperatureDirect() {}

  //-----------------------------------------------------------------------
  /// Return the temperature of the surface. This is different than the
  /// temperature near the surface which would be the lowest level of
  /// the temperature grid.
  //-----------------------------------------------------------------------

  virtual AutoDerivativeWithUnit<double>
  surface_temperature(int channel_index) const;
  virtual boost::shared_ptr<SurfaceTemperature> clone() const;

  virtual std::string sub_state_identifier() const
  { return "surface_temperature"; }

  std::string state_vector_name_i(int i) const;

private:
  Unit units;
  SurfaceTemperatureDirect() {}
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
typedef SubStateVectorArray<SurfaceTemperature> SubStateVectorArraySurfaceTemperature;
}

FP_EXPORT_KEY(SurfaceTemperatureDirect)
FP_EXPORT_KEY(SubStateVectorArraySurfaceTemperature);
#endif
