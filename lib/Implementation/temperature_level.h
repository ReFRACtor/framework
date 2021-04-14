#ifndef TEMPERATURE_LEVEL_H
#define TEMPERATURE_LEVEL_H

#include "temperature_imp_base.h"

namespace FullPhysics {
/****************************************************************//**
  This class maintains the temperature portion of the state. The
  temperature is retrieved through a set of levels related to pressures
  on the same grid.
*******************************************************************/

class TemperatureLevel: virtual public TemperatureImpBase {
public:
    TemperatureLevel(const blitz::Array<double, 1> Temp,
                     const boost::shared_ptr<Pressure>& Press,
                     boost::shared_ptr<StateMapping> Map = boost::make_shared<StateMappingLinear>());

    virtual ~TemperatureLevel() = default;

    virtual void print(std::ostream& Os) const;

    virtual boost::shared_ptr<Temperature> clone() const;

    virtual std::string sub_state_identifier() const
    {
        return "temperature_level";
    }

    virtual std::string state_vector_name_i(int i) const;

    //-----------------------------------------------------------------------
    /// Full temperature associate with the pressure profile, values are in 
    /// Kelvin
    //-----------------------------------------------------------------------

    virtual blitz::Array<double, 1> temperature_profile() const
    { return mapping->mapped_state(coeff).value(); }

    //-----------------------------------------------------------------------
    /// Pressure levels that serve as the grid for the temperature values in
    /// units of Pascals
    //-----------------------------------------------------------------------
    virtual blitz::Array<double, 1> pressure_profile() const
    { return press->pressure_grid().value.value(); }

    virtual ArrayWithUnit<double, 1> important_pressure_level() const
    {
        return ArrayWithUnit<double, 1>(pressure_profile(), units::Pa);
    }
protected:
    void calc_temperature_grid() const;
private:
    TemperatureLevel() = default;
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
};
}

FP_EXPORT_KEY(TemperatureLevel);

#endif
