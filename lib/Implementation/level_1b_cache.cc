#include "level_1b_cache.h"
#include "fp_serialize_support.h"
using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION

template<class Archive>
void Level1bCache::serialize(Archive & ar, const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Level1b)
    & FP_NVP(lat) & FP_NVP(lon) & FP_NVP(szen) & FP_NVP(sazm)
    & FP_NVP(solzen) & FP_NVP(solazm) & FP_NVP(alt) & FP_NVP(rvel)
    & FP_NVP(stk_coeff) & FP_NVP(samp_grid) & FP_NVP(tm)
    & FP_NVP(rad);
}

FP_IMPLEMENT(Level1bCache);
#endif

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
///
/// Note we leave off the units, since this gets called from python
/// where passing an array without units is easier. The angles are all
/// degrees. The altitude should be meters and relative velocity m/s.
/// Certainly the degrees coming in are unlikely to be anything other
/// than deg. I suppose the altitude and relative velocity could be in
/// other units - if this becomes an issue we can add in the units to
/// the argument list. We leave the units on for the Samp_grid and
/// Rad, these tend to vary so we really do want the units here.
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
    lat.push_back(DoubleWithUnit(Lat(i), "deg"));
    lon.push_back(DoubleWithUnit(Lon(i), "deg"));
    szen.push_back(DoubleWithUnit(Sounding_zenith(i), "deg"));
    sazm.push_back(DoubleWithUnit(Sounding_azimuth(i), "deg"));
    solzen.push_back(DoubleWithUnit(Solar_zenith(i), "deg"));
    solazm.push_back(DoubleWithUnit(Solar_azimuth(i), "deg"));
    alt.push_back(DoubleWithUnit(Altitude(i), "m"));
    rvel.push_back(DoubleWithUnit(Relative_velocity(i), "m/s"));
    blitz::Array<double, 1> stk(Stokes_coeff(i, ra));
    stk_coeff.push_back(stk.copy());
    samp_grid.push_back(Samp_grid[i]->clone());
    tm.push_back(Time(*Tm[i]));
    rad.push_back(Rad[i]->clone());
  }
}
