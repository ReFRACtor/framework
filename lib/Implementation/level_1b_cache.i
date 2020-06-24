// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "fp_common.i"
%{
#include "level_1b_cache.h"
%}
%base_import(level_1b)
%import "spectral_range.i"
%import "spectral_domain.i"
%fp_shared_ptr(FullPhysics::Level1bCache);

namespace FullPhysics {
class Level1bCache : public Level1b {
public:
  Level1bCache(const Level1b& L1_in);
  Level1bCache(const blitz::Array<double, 1>& Lat,
	       const blitz::Array<double, 1>& Lon,
	       const blitz::Array<double, 1>& Sounding_zenith,	       
	       const blitz::Array<double, 1>& Sounding_azimuth,
	       const blitz::Array<double, 1>& Solar_zenith,	       
	       const blitz::Array<double, 1>& Solar_azimuth,
	       const blitz::Array<double, 1>& Altitude,
	       const blitz::Array<double, 1>& Relative_velocity,
	       const blitz::Array<double, 2>& Stokes_coeff,
	       const std::vector<boost::shared_ptr<SpectralDomain> >& Samp_grid,
	       const std::vector<boost::shared_ptr<Time> >& Tm,
	       const std::vector<boost::shared_ptr<SpectralRange> >& Rad);
  %python_attribute(number_spectrometer, virtual int);
  virtual DoubleWithUnit latitude(int i) const;
  void set_latitude(int i, const DoubleWithUnit& V);
  virtual DoubleWithUnit longitude(int i) const;
  void set_longitude(int i, const DoubleWithUnit& V);
  virtual DoubleWithUnit sounding_zenith(int i) const;
  void set_sounding_zenith(int i, const DoubleWithUnit& V);
  virtual DoubleWithUnit sounding_azimuth(int i) const;
  void set_sounding_azimuth(int i, const DoubleWithUnit& V);
  virtual blitz::Array<double, 1> stokes_coefficient(int i) const;
  void set_stokes_coefficient(int i, const blitz::Array<double, 1>& V);
  virtual DoubleWithUnit solar_zenith(int i) const;
  void set_solar_zenith(int i, const DoubleWithUnit& V);
  virtual DoubleWithUnit solar_azimuth(int i) const;
  void set_solar_azimuth(int i, const DoubleWithUnit& V);
  virtual DoubleWithUnit altitude(int i) const;
  void set_altitude(int i, const DoubleWithUnit& V);
  virtual DoubleWithUnit relative_velocity(int i) const;
  void set_relative_velocity(int i, const DoubleWithUnit& V);  
  virtual SpectralDomain sample_grid(int i) const;
  void set_sample_grid(int i, const SpectralDomain& V);  
  virtual Time time(int i) const;
  void set_time(int i, const Time& V);
  virtual SpectralRange radiance(int Spec_index) const;
  void set_radiance(int i, const SpectralRange& V);
  void set_radiance(int i, const SpectralRange& V,
		    const std::vector<int>& Plist);
  %pickle_serialization();
};
}
