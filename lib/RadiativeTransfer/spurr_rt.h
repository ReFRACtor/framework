#ifndef SPURR_RT_H
#define SPURR_RT_H

#include "spurr_driver.h"
#include "radiative_transfer_single_wn.h"
#include "rt_atmosphere.h"
#include <boost/noncopyable.hpp>

namespace FullPhysics {

/****************************************************************//**
  Abstract Interface for Rt classes based on Spurr driver 
  implementations
 *******************************************************************/

class SpurrRt : public RadiativeTransferSingleWn,
                public Observer<RtAtmosphere>,
                public boost::noncopyable {
public:

  SpurrRt(const boost::shared_ptr<RtAtmosphere>& Atm,
          const boost::shared_ptr<StokesCoefficient>& Stokes_coef,
          const blitz::Array<double, 1>& Sza, 
          const blitz::Array<double, 1>& Zen, 
          const blitz::Array<double, 1>& Azm,
          bool do_solar = true,
          bool do_thermal = false);
  
  //-----------------------------------------------------------------------
  /// For performance, we cache some data as we calculate it. This
  /// becomes stale when the Atmosphere is changed, so we observe atm
  /// and mark the cache when it changes. 
  //-----------------------------------------------------------------------
  void notify_update(const RtAtmosphere& UNUSED(atm)) { alt_spec_index_cache = -1; }

  /// Number of stokes in returned stokes values
  /// Note that LIDORT will only ever calculate the first stoke index for I,
  virtual int number_stokes() const { return stokes_coef->stokes_coefficient().cols(); }

  /// Number of quadtature streams in the cosine half space
  virtual int number_stream() const = 0;

  /// Number of moments for scattering matrix.
  virtual int number_moment() const = 0;

  /// Integer representing the surface type using the LIDORT indexing nomenclature
  virtual int surface_type() const { return surface_type_int; }

  virtual void print(std::ostream& Os, bool Short_form = false) const;

  virtual blitz::Array<double, 1> stokes_single_wn(double Wn, int Spec_index, const boost::shared_ptr<OpticalProperties>& Opt_prop = NULL) const;
  virtual ArrayAd<double, 1> stokes_and_jacobian_single_wn(double Wn, int Spec_index, const boost::shared_ptr<OpticalProperties>& Opt_prop = NULL) const;

  //-----------------------------------------------------------------------
  /// The various SpurrRtDriver have a very large state. Most of this
  /// is just internal to the RT code (e.g, LIDORT), and we don't normally
  /// save all this. Instead, we save the little bit at the SpurrRt
  /// and derived class level, and then recreate the SpurrRtDriver
  /// when we load the serialization.
  ///
  /// If you are doing something special where you have directly
  /// modified the SpurrRtDriver (e.g., directly changed
  /// LidortBrdfDriver), then you can elect to save the complete state.
  ///
  /// Note that this is pretty sizable, something like 80MB, but if
  /// you need it this can be pretty useful.
  ///
  /// By default, we don't save the full state. Change this variable
  /// to "true" to save the complete state.
  //-----------------------------------------------------------------------

  static bool serialize_full_state;

  virtual std::vector<boost::shared_ptr<GenericObject> >
  subobject_list() const
  { std::vector<boost::shared_ptr<GenericObject> > res =
      RadiativeTransferSingleWn::subobject_list();
    res.push_back(rt_driver_);
    return res;
  }
protected:

  int surface_type_int;
  bool do_solar_sources, do_thermal_emission;

  blitz::Array<double, 1> sza, zen, azm;

  //-----------------------------------------------------------------------
  /// Note for serialization on load if this is null the derived class
  /// should create this to fill it in.
  //-----------------------------------------------------------------------

  boost::shared_ptr<SpurrRtDriver> rt_driver_;
  
  // Last index we updates the altitude/geometry for.
  mutable int alt_spec_index_cache, geo_spec_index_cache;
  virtual void update_altitude(int spec_index) const;
  virtual void update_geometry(int spec_index) const;
  SpurrRt() : alt_spec_index_cache(-1), geo_spec_index_cache(-1) {}
private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
  template<class Archive>
  void save(Archive & ar, const unsigned int version) const;
  template<class Archive>
  void load(Archive & ar, const unsigned int version);
};
}

FP_EXPORT_KEY(SpurrRt);
#endif

