#ifndef TEMPERATURE_FIXED_LEVEL_H
#define TEMPERATURE_FIXED_LEVEL_H
#include "pressure_level_input.h"
#include "temperature_imp_base.h"
#include "linear_interpolate.h"

namespace FullPhysics {
/****************************************************************//**
  This class maintains the temperature portion of the state. This
  particular implementation has a fixed set of pressure levels, with
  only the surface pressure changing.
*******************************************************************/
class TemperatureFixedLevel: public TemperatureImpBase {
public:
  TemperatureFixedLevel(const blitz::Array<bool, 1>& Flag_temp, 
	bool Flag_offset, const blitz::Array<double, 1>& Temp,
	double T_offset, 
        const boost::shared_ptr<Pressure>& Press,
        const boost::shared_ptr<PressureLevelInput>& Press_level);
  virtual ~TemperatureFixedLevel() {}
  ArrayAd<double, 1> temperature_levels() const;
  virtual void print(std::ostream& Os) const;

  virtual boost::shared_ptr<Temperature> clone() const;
  virtual std::string sub_state_identifier() const { return "temperature_levels"; }
  virtual std::string state_vector_name_i(int i) const;

//-----------------------------------------------------------------------
/// Temperature offset.
//-----------------------------------------------------------------------
  double temperature_offset() const {return coefficient()(0).value(); }
  double temperature_offset_uncertainty() const;
protected:
  void calc_temperature_grid() const;
private:
  boost::shared_ptr<PressureLevelInput> press_level;
  // Range of coefficient that contains temperature.
  blitz::Range temperature_range() const 
  {return blitz::Range(1, coefficient().rows() - 1);}
  /// Cache these variables, so the vmr function can access this data
  mutable std::vector<AutoDerivative<double> > tlist;
  mutable std::vector<AutoDerivative<double> > plist;
  typedef LinearInterpolate<AutoDerivative<double>, AutoDerivative<double> >
    lin_type;
  mutable boost::shared_ptr<lin_type> lin;
};
}
#endif
