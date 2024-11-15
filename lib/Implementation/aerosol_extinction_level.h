#ifndef AEROSOL_EXTINCTION_LEVEL_H
#define AEROSOL_EXTINCTION_LEVEL_H

#include "aerosol_extinction_imp_base.h"
#include "state_mapping_linear.h"
#include <boost/lexical_cast.hpp>

namespace FullPhysics {
/****************************************************************//**
  This class maps the state vector to the aerosol extinction on each
  level.

  This implementation just gets the extinction coefficient for each
  level from the state vector.
*******************************************************************/
class AerosolExtinctionLevel : virtual public AerosolExtinctionImpBase {
public:
//-----------------------------------------------------------------------
/// Constructor.
/// \param Press The pressure to use
/// \param Aext The aerosol extinction value.
/// \param Aerosol_name The name of the aerosol. This is used to
///   generate the state vector name metadata, so it should be
///   whatever is convenient.
/// \param in_map (optional) StateMapping that describes the relationship
///   between the forward model and retrieval views of Aerosol extinction.
///   Defaults to linear mapping.
//-----------------------------------------------------------------------

  AerosolExtinctionLevel(const boost::shared_ptr<Pressure>& Press,
                         const blitz::Array<double, 1>& Aext,
                         const std::string& Aerosol_name,
                         boost::shared_ptr<StateMapping> in_map = boost::make_shared<StateMappingLinear>())
  {
    init(Aerosol_name, Aext, Press, in_map);
  }
  virtual ~AerosolExtinctionLevel() {}

  virtual boost::shared_ptr<AerosolExtinction> clone() const;

  virtual std::string sub_state_identifier() const { return "aerosol_extinction/" + mapping->name() +
          "/" + aerosol_name(); }

  virtual std::string state_vector_name_i(int i) const
  { return "Aerosol " + aerosol_name() + " " + mapping->name() + " Aerosol Ext for Press Lvl " +
      boost::lexical_cast<std::string>(i + 1); }

  virtual std::string model_short_name() const { return "profile_" + mapping->name(); }
  virtual void print(std::ostream& Os) const;
protected:
  virtual void calc_aerosol_extinction() const;
  AerosolExtinctionLevel() {}
private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
}

FP_EXPORT_KEY(AerosolExtinctionLevel);
#endif
