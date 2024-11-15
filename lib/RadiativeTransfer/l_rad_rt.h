#ifndef L_RAD_RT_H
#define L_RAD_RT_H

#include "radiative_transfer_single_wn.h"
#include "rt_atmosphere.h"
#include "spectral_bound.h"
#include <boost/noncopyable.hpp>
#include <blitz/array.h>
#include "l_rad_driver.h"

namespace FullPhysics {
/****************************************************************//**
  Base class for LRadRt. This is basically everything that is not
  actually part of the RadiativeTransfer. We separate this out so
  we can have both LRadRt (a RadiativeTransferSingleWn) and 
  LRadRtFixedStokesCoefficient (a RadiativeTransferFixedStokesCoefficient).
*******************************************************************/

class LRadRtBase : public Observer<RtAtmosphere>,
		   public boost::noncopyable {
public:
  // These match numbers used internal to the Fortran code for LRad
  enum BrdfType {
		 LAMBERTIAN  = 1,
		 COXMUNK     = 2,
		 BPDFVEGN    = 3,
		 BPDFSOIL   = 4
  };
  
  LRadRtBase(const boost::shared_ptr<RtAtmosphere>& Atm,
	     const SpectralBound& Spec_bound,
	     const blitz::Array<double, 1>& Sza, 
	     const blitz::Array<double, 1>& Zen, 
	     const blitz::Array<double, 1>& Azm, 
	     bool Pure_nadir,
	     int Number_stokes,
	     bool Use_first_order_scatt_calc = true,
	     bool Do_second_order = false,
	     int Number_stream = 4,
	     double Spectrum_spacing = 0.01,
	     const LRadDriver::PsMode ps_mode = LRadDriver::DETECT,
	     bool Use_tms_correction = true);
  virtual ~LRadRtBase() {}
  virtual int number_stokes() const { return driver->number_stokes(); }
  virtual int number_stream() const { return driver->number_stream(); }

  /// Returns an integer with l_rad's representation of surface type
  virtual int surface_type() const;

  /// Return an interpolated z_matrix value for use in offline testing
  ArrayAd<double, 2> interp_z_matrix(double Wn) {
    return (*zmat_interpolate)(Wn);
  }

  //-----------------------------------------------------------------------
  /// For performance, we cache some data as we calculate it. This
  /// becomes stale when the Atmosphere is changed, so we observe atm
  /// and mark the cache when it changes. 
  //-----------------------------------------------------------------------
  void notify_update(const RtAtmosphere& UNUSED(atm))
  {
    alt_spec_index_cache = -1;
    driver->surface_type(surface_type());
  }

  blitz::Array<double, 1> stokes_corr
  (double Wn, int Spec_index, bool have_rt,
   const boost::shared_ptr<OpticalProperties>& Opt_prop = NULL) const;
  ArrayAd<double, 1> stokes_and_jacobian_corr
  (double Wn, int Spec_index, bool have_rt,
   const boost::shared_ptr<OpticalProperties>& Opt_prop = NULL) const;
  
  const boost::shared_ptr<LRadDriver>& l_rad_driver() const { return driver; };
  virtual void print(std::ostream& Os, bool Short_form = false) const;
protected:
  void initialize(const SpectralBound& Spec_bound, double Spectrum_spacing);

  /// Set up interpolated z-matrix giving it the ends of the bands, or window
  /// Using this reduces the need for costly calculations for each spectral point
  void setup_z_matrix_interpol(const double wmin, const ArrayAd<double, 3>& pf_min, const double wmax, const ArrayAd<double, 3>& pf_max) const;

  // Interpolation variables for z-matrix
  mutable boost::shared_ptr<LinearInterpolate2Point<double,
						    blitz::Array<double, 2> > > zmat_interpolate;
  mutable boost::shared_ptr<LinearInterpolate2Point<double,
						    blitz::Array<double, 3> > > l_zmat_interpolate;

  // Helper methods for stokes calculation methods
  ArrayAd<double, 2> get_z_matrix(const double Wn, int Spec_index, const boost::shared_ptr<OpticalProperties>& Opt_prop = NULL) const;
  void apply_jacobians(double Wn, int Spec_index, ArrayAd<double, 1>& stokes, const blitz::Array<double, 3>& jac_atm, const blitz::Array<double, 2>& jac_surf, const boost::shared_ptr<OpticalProperties>& Opt_prop = NULL) const;

  boost::shared_ptr<RtAtmosphere> atm_;
  boost::shared_ptr<Ground> ground_;
  
  // Control how to use the l_rad_driver
  bool use_first_order_scatt_calc;
  bool do_second_order;

  // Scene geometry
  blitz::Array<double, 1> sza, zen, azm;

  // Mininum and maximum wavenumbers per spectrometer for z-matrix interpolation
  std::vector<double> wmin, wmax;

  // Driver for l_rad Fortran code
  boost::shared_ptr<LRadDriver> driver;

  // Last index we updates the altitude for.
  mutable int alt_spec_index_cache;
  void update_altitude(int spec_index) const;
  LRadRtBase() : alt_spec_index_cache(-1) {}
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
  
/****************************************************************//**
  This class drives the LRAD code, which gives a polarization
  correction to scalar intensity and jacobians. This can also be used
  on its own to provide a single scatter approximation to the
  RadiativeTransfer (i.e., without also running LIDORT).
*******************************************************************/

class LRadRt : public RadiativeTransferSingleWn,
	       public LRadRtBase {
public:
  LRadRt(const boost::shared_ptr<RadiativeTransferSingleWn>& Rt,
	 const SpectralBound& Spec_bound,
	 const blitz::Array<double, 1>& Sza, 
	 const blitz::Array<double, 1>& Zen, 
	 const blitz::Array<double, 1>& Azm, 
	 bool Pure_nadir,
	 bool Use_first_order_scatt_calc = true,
	 bool Do_second_order = false,
	 double Spectrum_spacing = 0.01,
	 const LRadDriver::PsMode ps_mode = LRadDriver::DETECT);
  LRadRt(const boost::shared_ptr<StokesCoefficient>& Stokes_coef,
	 const boost::shared_ptr<RtAtmosphere>& Atm,
	 const SpectralBound& Spec_bound,
	 const blitz::Array<double, 1>& Sza, 
	 const blitz::Array<double, 1>& Zen, 
	 const blitz::Array<double, 1>& Azm, 
	 bool Pure_nadir,
	 int Number_stokes,
	 bool Do_second_order = false,
	 int Number_stream = 4,
	 double Spectrum_spacing = 0.01,
	 const LRadDriver::PsMode ps_mode = LRadDriver::DETECT);
  virtual int number_stokes() const { return LRadRtBase::number_stokes(); }
  virtual int number_stream() const { return LRadRtBase::number_stream(); }
  virtual blitz::Array<double, 1> stokes_single_wn
  (double Wn, int Spec_index,
   const boost::shared_ptr<OpticalProperties>& Opt_prop = NULL) const;
  virtual ArrayAd<double, 1> stokes_and_jacobian_single_wn
  (double Wn, int Spec_index,
   const boost::shared_ptr<OpticalProperties>& Opt_prop = NULL) const;
  const boost::shared_ptr<RadiativeTransferSingleWn>& radiative_transfer()
    const { return rt; }
  virtual void print(std::ostream& Os, bool Short_form = false) const;
  virtual std::vector<boost::shared_ptr<GenericObject> >
  subobject_list() const
  { std::vector<boost::shared_ptr<GenericObject> > res =
      RadiativeTransferSingleWn::subobject_list();
    res.push_back(rt);
    res.push_back(driver);
    return res;
  }
private:
  // Underlying RT that we are applying correction for.
  boost::shared_ptr<RadiativeTransferSingleWn> rt;
  LRadRt() {}
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};

/****************************************************************//**
  This is a variation of LRadRt that works with 
  RadiativeTransferFixedStokesCoefficient (e.g., PCARt).
*******************************************************************/

class LRadRtFixedStokesCoefficient :
    public RadiativeTransferFixedStokesCoefficient,
    public LRadRtBase {
public:
  LRadRtFixedStokesCoefficient
  (const boost::shared_ptr<RadiativeTransferFixedStokesCoefficient>& Rt,
   const boost::shared_ptr<RtAtmosphere>& Atm,
   int Number_stream,
   const SpectralBound& Spec_bound,
   const blitz::Array<double, 1>& Sza, 
   const blitz::Array<double, 1>& Zen, 
   const blitz::Array<double, 1>& Azm,
   bool Pure_nadir,
   bool Use_first_order_scatt_calc = true,
   bool Do_second_order = false,
   double Spectrum_spacing = 0.01,
   const LRadDriver::PsMode ps_mode = LRadDriver::DETECT);
  virtual int number_stokes() const { return LRadRtBase::number_stokes(); }
  virtual int number_stream() const { return LRadRtBase::number_stream(); }
  virtual blitz::Array<double, 2> stokes(const SpectralDomain& Spec_domain,
					 int Spec_index) const;
  virtual ArrayAd<double, 2> stokes_and_jacobian
  (const SpectralDomain& Spec_domain, int Spec_index) const;
  const boost::shared_ptr<RadiativeTransferFixedStokesCoefficient>&
  radiative_transfer() const { return rt; }
  virtual void print(std::ostream& Os, bool Short_form = false) const;
  virtual std::vector<boost::shared_ptr<GenericObject> >
  subobject_list() const
  { std::vector<boost::shared_ptr<GenericObject> > res =
      RadiativeTransferFixedStokesCoefficient::subobject_list();
    res.push_back(rt);
    res.push_back(driver);
    return res;
  }
private:
  // Underlying RT that we are applying correction for.
  boost::shared_ptr<RadiativeTransferFixedStokesCoefficient> rt;
  LRadRtFixedStokesCoefficient() {}
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
}

FP_EXPORT_KEY(LRadRtBase);
FP_EXPORT_KEY(LRadRt);
FP_EXPORT_KEY(LRadRtFixedStokesCoefficient);
FP_CLASS_VERSION(LRadRt, 1);
#endif
