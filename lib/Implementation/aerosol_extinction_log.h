#ifndef AEROSOL_EXTINCTION_LOG_H
#define AEROSOL_EXTINCTION_LOG_H

#include "aerosol_extinction_level.h"
#include "pressure.h"
#include <boost/lexical_cast.hpp>

namespace FullPhysics {
/****************************************************************//**
  This class maps the state vector to the aerosol extinction on each
  level.

  This implementation just gets the log of the extinction coefficient
  for each level from the state vector.
*******************************************************************/
class AerosolExtinctionLog : virtual public AerosolExtinctionLevel {
public:
//-----------------------------------------------------------------------
/// Constructor.
/// \param Press The pressure to use
/// \param Aext The aerosol extinction value.
/// \param Aerosol_name The name of the aerosol. This is used to
///   generate the state vector name metadata, so it should be
///   whatever is convenient.
//-----------------------------------------------------------------------

  AerosolExtinctionLog(const boost::shared_ptr<Pressure>& Press,
			  const blitz::Array<double, 1>& Aext,
			  const std::string& Aerosol_name);

  virtual ~AerosolExtinctionLog() = default;

  virtual boost::shared_ptr<AerosolExtinction> clone() const;
protected:
  AerosolExtinctionLog() {}
private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
}

FP_EXPORT_KEY(AerosolExtinctionLog);
#endif
