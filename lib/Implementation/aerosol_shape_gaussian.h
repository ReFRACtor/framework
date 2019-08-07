#ifndef AEROSOL_SHAPE_GAUSSIAN_H
#define AEROSOL_SHAPE_GAUSSIAN_H

#include "aerosol_extinction_level.h"
#include "pressure.h"
#include <boost/lexical_cast.hpp>

namespace FullPhysics {
/****************************************************************//**
  This class maps the state vector to aerosol extinction defined
  by a Gaussian parameterization.
*******************************************************************/
class AerosolShapeGaussian : public AerosolExtinctionLevel {
public:
//-----------------------------------------------------------------------
/// Constructor.
/// \param Press Pressure object used for defining pressure levels
/// \param Flag Boolean flag indicating which coefficients are to be set by
///   the state vector. A value of false means it is held fixed
///   when the state vector changes.
/// \param Coeffs The total aerosol optical depth, center pressure
///   and pressure width coefficients. The total optical depth is
///   defined in linear or log terms. The pressure coefficents are relative
///   to the surface pressure and are hence dimensionless.
/// \param Aerosol_name The name of the aerosol. This is used to
///   generate the state vector name metadata, so it should be
///   whatever is convenient.
/// \param Linear_AOD Specifies whether the total aod cofficient
///   is in linear space if true, logarithmic space if false.
//-----------------------------------------------------------------------

AerosolShapeGaussian(const boost::shared_ptr<Pressure>& Press,
		     const blitz::Array<bool, 1>& Flag, 
		     const blitz::Array<double, 1>& Coeffs,
		     const std::string& Aerosol_name,
		     const bool Linear_AOD);

virtual boost::shared_ptr<AerosolExtinction> clone() const
  { return clone(press->clone()); }
  virtual boost::shared_ptr<AerosolExtinction> clone
(const boost::shared_ptr<Pressure>& P) const;
  virtual ~AerosolShapeGaussian() {}
private:
  static const double min_aod;
};
}
#endif
