// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "fp_common.i"
%{
#include "aerosol_property_rh_hdf.h"
#include "absorber.h"
#include "temperature.h"
#include "altitude.h"
%}
%base_import(aerosol_property_imp_base)
%import "hdf_file.i"
%import "relative_humidity.i"
%fp_shared_ptr(FullPhysics::AerosolPropertyRhHdf)
namespace FullPhysics {
class AerosolPropertyRhHdf : public AerosolPropertyImpBase {
public:
  AerosolPropertyRhHdf(const HdfFile& F, const std::string& Group_name,
		       const boost::shared_ptr<Pressure>& Press,
		       const boost::shared_ptr<RelativeHumidity>& Rh);
  virtual boost::shared_ptr<AerosolProperty> clone() const;
  virtual boost::shared_ptr<AerosolProperty> 
  clone(const boost::shared_ptr<Pressure>& Press,
	const boost::shared_ptr<RelativeHumidity>& Rh) const;
  virtual ArrayAd<double, 1> extinction_coefficient_each_layer(double wn) const;
  virtual ArrayAd<double, 1> scattering_coefficient_each_layer(double wn) const;
  virtual ArrayAd<double, 3> 
  phase_function_moment_each_layer(double wn, int nmom = -1, 
				   int nscatt = -1) const;

};
}

