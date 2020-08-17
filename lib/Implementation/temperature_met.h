#ifndef TEMPERATURE_MET_H
#define TEMPERATURE_MET_H
#include "temperature_offset.h"
#include "meteorology.h"

namespace FullPhysics {
/****************************************************************//**
  This class maintains the temperature portion of the state. This
  particular implementation uses the temperature from MET file
  (interpolated to the current pressure grid), along with an offset.
*******************************************************************/
class TemperatureMet: virtual public TemperatureOffset {
public:
  TemperatureMet(const boost::shared_ptr<Meteorology>& Met_file,
                 const boost::shared_ptr<Pressure>& Press,
                 double Temp_offset);
  virtual ~TemperatureMet() {}
  virtual void print(std::ostream& Os) const;

  virtual boost::shared_ptr<Temperature> clone() const;

//-----------------------------------------------------------------------
/// Temperature from MET, used to write to output file
//-----------------------------------------------------------------------

  virtual blitz::Array<double, 1> temperature_profile() const
      { return met->temperature(); }

//-----------------------------------------------------------------------
/// Pressure levels that temperature is on from MET, used to write
/// to output file 
//-----------------------------------------------------------------------

  virtual blitz::Array<double, 1> pressure_profile() const
      { return met->pressure_levels(); } 

private:
  boost::shared_ptr<Meteorology> met;
  TemperatureMet() {}
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
}
FP_EXPORT_KEY(TemperatureMet);
#endif
