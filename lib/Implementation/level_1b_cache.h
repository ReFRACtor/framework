#ifndef LEVEL_1B_CACHE_H
#define LEVEL_1B_CACHE_H
#include "level_1b.h"
#include <vector>
namespace FullPhysics {
/****************************************************************//**
  This is a Level1b implementation that just saves the data read from
  another Level1b object, or the data can be directly passed in.
  We then allow these values to be changed if
  desired. This can be useful when setting up special run in Python, 
  among other uses.
*******************************************************************/
class Level1bCache: public Level1b {
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
  virtual ~Level1bCache() {}
  virtual void print(std::ostream& Os) const {Os << "Level1bCache";}
  virtual int number_spectrometer() const { return (int) lat.size(); }
  virtual DoubleWithUnit latitude(int i) const 
  { 
    range_check(i, 0, number_spectrometer());
    return lat[i];
  }

//-----------------------------------------------------------------------
/// Change value.
//-----------------------------------------------------------------------

  void set_latitude(int i, const DoubleWithUnit& V)
  { 
    range_check(i, 0, number_spectrometer());
    lat[i] = V;
  }
  virtual DoubleWithUnit longitude(int i) const
  { 
    range_check(i, 0, number_spectrometer());
    return lon[i];
  }

//-----------------------------------------------------------------------
/// Change value.
//-----------------------------------------------------------------------

  void set_longitude(int i, const DoubleWithUnit& V)
  { 
    range_check(i, 0, number_spectrometer());
    lon[i] = V;
  }
  virtual DoubleWithUnit sounding_zenith(int i) const
  { 
    range_check(i, 0, number_spectrometer());
    return szen[i];
  }

//-----------------------------------------------------------------------
/// Change value.
//-----------------------------------------------------------------------

  void set_sounding_zenith(int i, const DoubleWithUnit& V)
  { 
    range_check(i, 0, number_spectrometer());
    szen[i] = V;
  }
  virtual DoubleWithUnit sounding_azimuth(int i) const
  { 
    range_check(i, 0, number_spectrometer());
    return sazm[i];
  }

//-----------------------------------------------------------------------
/// Change value.
//-----------------------------------------------------------------------

  void set_sounding_azimuth(int i, const DoubleWithUnit& V)
  { 
    range_check(i, 0, number_spectrometer());
    sazm[i] = V;
  }
  virtual DoubleWithUnit solar_zenith(int i) const
  { 
    range_check(i, 0, number_spectrometer());
    return solzen[i];
  }

//-----------------------------------------------------------------------
/// Change value.
//-----------------------------------------------------------------------

  void set_solar_zenith(int i, const DoubleWithUnit& V)
  { 
    range_check(i, 0, number_spectrometer());
    solzen[i] = V;
  }
  virtual DoubleWithUnit solar_azimuth(int i) const
  { 
    range_check(i, 0, number_spectrometer());
    return solazm[i];
  }

//-----------------------------------------------------------------------
/// Change value.
//-----------------------------------------------------------------------

  void set_solar_azimuth(int i, const DoubleWithUnit& V)
  { 
    range_check(i, 0, number_spectrometer());
    solazm[i] = V;
  }
  virtual DoubleWithUnit altitude(int i) const
  { 
    range_check(i, 0, number_spectrometer());
    return alt[i];
  }

//-----------------------------------------------------------------------
/// Change value.
//-----------------------------------------------------------------------

  void set_altitude(int i, const DoubleWithUnit& V) 
  { 
    range_check(i, 0, number_spectrometer());
    alt[i] = V;
  }
  virtual DoubleWithUnit relative_velocity(int i) const
  { 
    range_check(i, 0, number_spectrometer());
    return rvel[i];
  }

//-----------------------------------------------------------------------
/// Change value.
//-----------------------------------------------------------------------

  void set_relative_velocity(int i, const DoubleWithUnit& V)
  { 
    range_check(i, 0, number_spectrometer());
    rvel[i] = V;
  }
  virtual blitz::Array<double, 1> stokes_coefficient(int i) const
  {
    range_check(i, 0, number_spectrometer());
    return stk_coeff[i];
  }    

//-----------------------------------------------------------------------
/// Change value.
//-----------------------------------------------------------------------

  void set_stokes_coefficient(int i, const blitz::Array<double, 1>& V)
  { 
    range_check(i, 0, number_spectrometer());
    stk_coeff[i].reference(V.copy());
  }

  SpectralDomain sample_grid(int i) const
  {
    range_check(i, 0, number_spectrometer());
    return samp_grid[i];
  }

//-----------------------------------------------------------------------
/// Change value.
//-----------------------------------------------------------------------

  void set_sample_grid(int i, const SpectralDomain& V)
  {
    range_check(i, 0, number_spectrometer());
    samp_grid[i] = V;
  }

  virtual Time time(int i) const
  {
    range_check(i, 0, number_spectrometer());
    return tm[i];
  }

//-----------------------------------------------------------------------
/// Change value.
//-----------------------------------------------------------------------

  void set_time(int i, const Time& V)
  { 
    range_check(i, 0, number_spectrometer());
    tm[i] = V;
  }
  virtual SpectralRange radiance(int i) const
  { 
    range_check(i, 0, number_spectrometer());
    return rad[i];
  }

//-----------------------------------------------------------------------
/// Change value.
//-----------------------------------------------------------------------

  void set_radiance(int i, const SpectralRange& V)
  { 
    range_check(i, 0, number_spectrometer());
    rad[i] = V;
  }

//-----------------------------------------------------------------------
/// Change value, but only for a subset of pixels. This might come
/// from the ForwardModelSpectralGrid for example.
//-----------------------------------------------------------------------

 void set_radiance(int i, const SpectralRange& V,
		    const std::vector<int>& Plist)
  { 
    range_check(i, 0, number_spectrometer());
    if(V.data().rows() != (int) Plist.size())
      throw Exception("Spectral Range is not the same size as Plist");
    blitz::Array<double, 1> rnew = rad[i].data().copy();
    blitz::Array<double, 1> unew = rad[i].uncertainty().copy();
    for(int j = 0; j < (int) Plist.size(); ++j) {
      rnew(Plist[j]) = V.data()(j);
      unew(Plist[j]) = V.uncertainty()(j);
    }
    rad[i] = SpectralRange(rnew, V.units(), unew);
  }

private:
  std::vector<DoubleWithUnit> lat, lon, szen, sazm, solzen, solazm, alt,
						    rvel;
  std::vector<blitz::Array<double, 1> > stk_coeff;
  std::vector<SpectralDomain> samp_grid;
  std::vector<Time> tm;
  std::vector<SpectralRange> rad;
  Level1bCache() {}
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
}

FP_EXPORT_KEY(Level1bCache);
#endif
