#include "level_1b_cache.h"
using namespace FullPhysics;

//-----------------------------------------------------------------------
/// Constructor.
//-----------------------------------------------------------------------

Level1bCache::Level1bCache(const Level1b& L1_in)
{
  for(int i = 0; i < L1_in.number_spectrometer(); ++i) {
    lat.push_back(L1_in.latitude(i));
    lon.push_back(L1_in.longitude(i));
    szen.push_back(L1_in.sounding_zenith(i));
    sazm.push_back(L1_in.sounding_azimuth(i));
    solzen.push_back(L1_in.solar_zenith(i));
    solazm.push_back(L1_in.solar_azimuth(i));
    alt.push_back(L1_in.altitude(i));
    rvel.push_back(L1_in.relative_velocity(i));
    stk_coeff.push_back(L1_in.stokes_coefficient(i).copy());
    samp_grid.push_back(L1_in.sample_grid(i));
    tm.push_back(L1_in.time(i));
    rad.push_back(L1_in.radiance(i).clone());
  }
}

//-----------------------------------------------------------------------
/// Constructor taking data.
//-----------------------------------------------------------------------

Level1bCache::Level1bCache
(const blitz::Array<double, 1>& Lat,
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
 const std::vector<boost::shared_ptr<SpectralRange> >& Rad)
{
  blitz::Range ra = blitz::Range::all();
  if(Lat.rows() != Lon.rows() ||
     Lat.rows() != Sounding_zenith.rows() ||
     Lat.rows() != Sounding_azimuth.rows() ||
     Lat.rows() != Solar_zenith.rows() ||
     Lat.rows() != Solar_azimuth.rows() ||
     Lat.rows() != Altitude.rows() ||
     Lat.rows() != Relative_velocity.rows() ||
     Lat.rows() != Stokes_coeff.rows() ||
     Lat.rows() != (int) Samp_grid.size() ||
     Lat.rows() != (int) Tm.size() ||
     Lat.rows() != (int) Rad.size())
    throw Exception("All the arguments need to be the same size");

  for(int i = 0; i < Lat.rows(); ++i) {
    lat.push_back(Lat(i));
    lon.push_back(Lon(i));
    szen.push_back(Sounding_zenith(i));
    sazm.push_back(Sounding_azimuth(i));
    solzen.push_back(Solar_zenith(i));
    solazm.push_back(Solar_azimuth(i));
    alt.push_back(Altitude(i));
    rvel.push_back(Relative_velocity(i));
    blitz::Array<double, 1> stk(Stokes_coeff(i, ra));
    stk_coeff.push_back(stk.copy());
    samp_grid.push_back(Samp_grid[i]->clone());
    tm.push_back(Time(*Tm[i]));
    rad.push_back(Rad[i]->clone());
  }
}
