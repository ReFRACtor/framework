#ifndef AEROSOL_EXTINCTION_LINEAR_H
#define AEROSOL_EXTINCTION_LINEAR_H

#include "aerosol_extinction_level.h"
#include <boost/lexical_cast.hpp>

namespace FullPhysics {
/****************************************************************//**
  This class maps the state vector to the aerosol extinction on each
  level.

  This implementation just gets the extinction coefficient for each
  level from the state vector.
*******************************************************************/
class AerosolExtinctionLinear : public AerosolExtinctionLevel {
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
//-----------------------------------------------------------------------

  AerosolExtinctionLinear(const boost::shared_ptr<Pressure>& Press,
			  const blitz::Array<bool, 1>& Flag, 
			  const blitz::Array<double, 1>& Aext,
			  const std::string& Aerosol_name);

  virtual ~AerosolExtinctionLinear() = default;
    
  virtual boost::shared_ptr<AerosolExtinction> clone() const;
};
}
#endif
