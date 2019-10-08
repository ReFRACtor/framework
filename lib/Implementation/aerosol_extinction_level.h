#ifndef AEROSOL_EXTINCTION_LEVEL_H
#define AEROSOL_EXTINCTION_LEVEL_H

#include "aerosol_extinction_imp_base.h"
#include "mapping.h"
#include "mapping_linear.h"
#include <boost/lexical_cast.hpp>

namespace FullPhysics {
/****************************************************************//**
  This class maps the state vector to the aerosol extinction on each
  level.

  This implementation just gets the extinction coefficient for each
  level from the state vector.
*******************************************************************/
class AerosolExtinctionLevel : public AerosolExtinctionImpBase {
public:
//-----------------------------------------------------------------------
/// Constructor.
/// \param Press The pressure to use
/// \param Flag Boolean flag indicating which levels are to be set by
///   the state vector. A value of false means the level is held fixed
///   when the state vector changes.
/// \param Aext The aerosol extinction value.
/// \param Aerosol_name The name of the aerosol. This is used to
///   generate the state vector name metadata, so it should be
///   whatever is convenient.
/// \param in_map (optional) Mapping that describes the relationship
///   between the forward model and retrieval views of Aerosol extinction.
///   Defaults to linear mapping.
//-----------------------------------------------------------------------

  AerosolExtinctionLevel(const boost::shared_ptr<Pressure>& Press,
			  const blitz::Array<bool, 1>& Flag, 
			  const blitz::Array<double, 1>& Aext,
			  const std::string& Aerosol_name,
			  boost::shared_ptr<Mapping> in_map = boost::make_shared<MappingLinear>())
  {
    bool Mark_according_to_press = false;
    int Pdep_start = 0;
    init(Aerosol_name, Aext, Flag, Press, Mark_according_to_press,
                               Pdep_start, in_map);
  }
  virtual ~AerosolExtinctionLevel() {}

  virtual boost::shared_ptr<AerosolExtinction> clone() const { return clone(press->clone()); }

  virtual boost::shared_ptr<AerosolExtinction> clone(const boost::shared_ptr<Pressure>& P) const;

  virtual std::string sub_state_identifier() const { return "aerosol_extinction/" + mapping->name() +
          "/" + aerosol_name(); }

  virtual std::string state_vector_name_i(int i) const
  { return "Aerosol " + aerosol_name() + " " + mapping->name() + " Aerosol Ext for Press Lvl " +
      boost::lexical_cast<std::string>(i + 1); }

  virtual std::string model_short_name() const { return "profile_" + mapping->name(); }
  virtual void print(std::ostream& Os) const;
protected:
  virtual void calc_aerosol_extinction() const;
};
}
#endif
