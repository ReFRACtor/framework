#ifndef TEMPERATURE_LEVEL_OFFSET_H
#define TEMPERATURE_LEVEL_OFFSET_H
#include "temperature_offset.h"

namespace FullPhysics {
/****************************************************************//**
  This class maintains the temperature portion of the state. This
  particular implementation uses the temperature from LEVEL file
  (interpolated to the current pressure grid), along with an offset.
*******************************************************************/
class TemperatureLevelOffset: virtual public TemperatureOffset {
public:
  TemperatureLevelOffset(const boost::shared_ptr<Pressure>& Press,
                         const blitz::Array<double, 1>& Temp_levels,
                         double Temp_offset);
  virtual ~TemperatureLevelOffset() {}
  virtual void print(std::ostream& Os) const;

  virtual boost::shared_ptr<Temperature> clone() const;

  //-----------------------------------------------------------------------
  /// Temperature from profile, used to write to output file
  //-----------------------------------------------------------------------
  blitz::Array<double, 1> temperature_profile() const { return temp_levels; }
  blitz::Array<double, 1> pressure_profile() const { return press->pressure_grid().value.value(); }

private:
  blitz::Array<double, 1> temp_levels;
  TemperatureLevelOffset() {}
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
}
FP_EXPORT_KEY(TemperatureLevelOffset);
#endif
