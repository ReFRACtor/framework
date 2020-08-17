#ifndef TEMPERATURE_OFFSET_H
#define TEMPERATURE_OFFSET_H
#include "temperature_imp_base.h"

namespace FullPhysics {
/****************************************************************//**
  This class maintains the temperature portion of the state. The
  temperature is retrieved through an offset to an initial guess
  profile.
*******************************************************************/
class TemperatureOffset: virtual public TemperatureImpBase {
public:
  TemperatureOffset(const boost::shared_ptr<Pressure>& Press,
                    double Temp_offset);
  virtual ~TemperatureOffset() {}
  virtual void print(std::ostream& Os) const;

  virtual boost::shared_ptr<Temperature> clone() const = 0;

  virtual std::string sub_state_identifier() const { return "temperature_offset"; }
  virtual std::string state_vector_name_i(int i) const;

  //-----------------------------------------------------------------------
  /// Temperature offset.
  //-----------------------------------------------------------------------
  double temperature_offset() const {return coefficient()(0).value(); }
  double temperature_offset_uncertainty() const;

  //-----------------------------------------------------------------------
  /// Temperature from ECMWF, used to write to output file
  //-----------------------------------------------------------------------

  virtual blitz::Array<double, 1> temperature_profile() const = 0;

  //-----------------------------------------------------------------------
  /// Pressure levels that temperature is on from ECMWF, used to write
  /// to output file 
  //-----------------------------------------------------------------------
  virtual blitz::Array<double, 1> pressure_profile() const = 0;

  virtual ArrayWithUnit<double, 1> important_pressure_level() const
  {
    return ArrayWithUnit<double, 1>(pressure_profile(), units::Pa);
  }
protected:
  void calc_temperature_grid() const;
  TemperatureOffset() {}
private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
}

FP_EXPORT_KEY(TemperatureOffset);

#endif
