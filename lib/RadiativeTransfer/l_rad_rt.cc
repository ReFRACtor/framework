#include "l_rad_rt.h"
#include "fp_serialize_support.h"
#include "ostream_pad.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void LRadRtBase::serialize(Archive & ar,
			const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ObserverRtAtmosphere)
    & FP_NVP_(atm)
    & FP_NVP_(ground)
    & FP_NVP(use_first_order_scatt_calc)
    & FP_NVP(do_second_order) & FP_NVP(sza) & FP_NVP(zen)
    & FP_NVP(azm) & FP_NVP(wmin) & FP_NVP(wmax)
    & FP_NVP(driver);
}

template<class Archive>
void LRadRt::serialize(Archive & ar,
			const unsigned int version)
{
  if(version == 0) {
    // We didn't have LRadRtBase originally, so use the older
    // interface for older data
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(RadiativeTransferSingleWn)
      & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ObserverRtAtmosphere)
      & FP_NVP_(ground)
      & FP_NVP(use_first_order_scatt_calc)
      & FP_NVP(do_second_order) & FP_NVP(sza) & FP_NVP(zen)
      & FP_NVP(azm) & FP_NVP(wmin) & FP_NVP(wmax) & FP_NVP(rt)
      & FP_NVP(driver);
    // Because of the inheritance structure, LRadRtBase doesn't have
    // access to the atm in the RadiativeTransferSingleWn. So we have
    // a copy, which we need to make explicit because the Old
    // serialization version didn't have the l_rad base.
    atm_ = atm;
    return;
  }
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(LRadRtBase)
    & BOOST_SERIALIZATION_BASE_OBJECT_NVP(RadiativeTransferSingleWn)
    & FP_NVP(rt);
}

template<class Archive>
void LRadRtFixedStokesCoefficient::serialize(Archive & ar,
		     const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(LRadRtBase)
    & BOOST_SERIALIZATION_BASE_OBJECT_NVP(RadiativeTransferFixedStokesCoefficient)
    & FP_NVP(rt);
}

FP_IMPLEMENT(LRadRtBase);
FP_IMPLEMENT(LRadRt);
FP_IMPLEMENT(LRadRtFixedStokesCoefficient);
#endif

#ifdef HAVE_LUA
#include "register_lua.h"
boost::shared_ptr<RadiativeTransfer>
l_rad_rt_create(const boost::shared_ptr<RadiativeTransfer>& Rt,
                const SpectralBound& Spec_bound,
                const blitz::Array<double, 1>& Sza,
                const blitz::Array<double, 1>& Zen,
                const blitz::Array<double, 1>& Azm,
                bool Pure_nadir,
                bool Use_first_order_scatt_calc,
                bool Do_second_order)
{
    boost::shared_ptr<RadiativeTransferSingleWn> rts =
        boost::dynamic_pointer_cast<RadiativeTransferSingleWn>(Rt);

    return boost::make_shared<LRadRt>(rts, Spec_bound, Sza, Zen, Azm,
				      Pure_nadir, Use_first_order_scatt_calc,
				      Do_second_order);
}

REGISTER_LUA_DERIVED_CLASS(LRadRt, RadiativeTransfer)
.scope
[
    luabind::def("create", &l_rad_rt_create)
]
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Constructor. Note that we can run in two modes, either we are
/// doing a polarization correction to an underlying multi-scatter RT
/// code, or we are just doing a single-scatter calculation in the
/// LRad alone. This constructor sets up for the mult-scatter
/// correction.
///
/// \param Rt RT to apply polarization correction to.
/// \param Spec_bound Spectral window bounds for each spectrometer.
/// \param Sza Solar zenith angle. This is in degrees, and should be
///     in the range 0 to 90, and have size number_spectrometer()
/// \param Zen Zenith angle (degrees), in range 0 to 90, and have size
///      number_spectrometer()
/// \param Azm Azimuth angle (degrees), in range 0 to 360, and have size
///      number_spectrometer()
/// \param Pure_nadir If true then scene is purely nadir
/// \param Use_first_order_scatt_calc Use the first order of scattering
///      calculation from LRad since it is a faster calculation than
///      by most other radiative transfers. If using two RTs together
///      the other RT must only calculate the multiple scattering component
///      of the stokes vector.
/// \param Do_second_order If true, we do second order corrections
/// \param Spectrum_spacing The spectrum spacing
/// \param ps_mode Which pseudo spherical mode to use: 
///        REGULAR, ENHANCED, PLANE_PARALLEL, DETECT
//-----------------------------------------------------------------------

LRadRt::LRadRt(const boost::shared_ptr<RadiativeTransferSingleWn>& Rt,
               const SpectralBound& Spec_bound,
               const blitz::Array<double, 1>& Sza,
               const blitz::Array<double, 1>& Zen,
               const blitz::Array<double, 1>& Azm,
               bool Pure_nadir,
               bool Use_first_order_scatt_calc,
               bool Do_second_order,
               double Spectrum_spacing,
               const LRadDriver::PsMode Ps_mode)
: RadiativeTransferSingleWn(Rt->stokes_coefficient(),
			    Rt->atmosphere_ptr()),
  LRadRtBase(Rt->atmosphere_ptr(), Spec_bound,
	     Sza, Zen, Azm, Pure_nadir, Rt->number_stokes(),
	     Use_first_order_scatt_calc,
	     Do_second_order, Rt->number_stream(),
	     Spectrum_spacing, Ps_mode,
	     true // Use TMS correction, have RT
	     ),
  rt(Rt)
{
}

//-----------------------------------------------------------------------
/// Constructor. Note that we can run in two modes, either we are
/// doing a polarization correction to an underlying multi-scatter RT
/// code, or we are just doing a single-scatter calculation in the
/// LRad alone. This constructor sets up for the mult-scatter
/// correction.
///
/// \param Rt RT to apply polarization correction to.
/// \param Atm Atmosphere
/// \param Spec_bound Spectral window bounds for each spectrometer.
/// \param Sza Solar zenith angle. This is in degrees, and should be
///     in the range 0 to 90, and have size number_spectrometer()
/// \param Zen Zenith angle (degrees), in range 0 to 90, and have size
///      number_spectrometer()
/// \param Azm Azimuth angle (degrees), in range 0 to 360, and have size
///      number_spectrometer()
/// \param Pure_nadir If true then scene is purely nadir
/// \param Use_first_order_scatt_calc Use the first order of scattering
///      calculation from LRad since it is a faster calculation than
///      by most other radiative transfers. If using two RTs together
///      the other RT must only calculate the multiple scattering component
///      of the stokes vector.
/// \param Do_second_order If true, we do second order corrections
/// \param Spectrum_spacing The spectrum spacing
/// \param ps_mode Which pseudo spherical mode to use: 
///        REGULAR, ENHANCED, PLANE_PARALLEL, DETECT
//-----------------------------------------------------------------------

LRadRtFixedStokesCoefficient::LRadRtFixedStokesCoefficient
(const boost::shared_ptr<RadiativeTransferFixedStokesCoefficient>& Rt,
 const boost::shared_ptr<RtAtmosphere>& Atm,
 int Number_stream,
 const SpectralBound& Spec_bound,
 const blitz::Array<double, 1>& Sza,
 const blitz::Array<double, 1>& Zen,
 const blitz::Array<double, 1>& Azm,
 bool Pure_nadir,
 bool Use_first_order_scatt_calc,
 bool Do_second_order,
 double Spectrum_spacing,
 const LRadDriver::PsMode Ps_mode)
: RadiativeTransferFixedStokesCoefficient(Rt->stokes_coefficient()),
  LRadRtBase(Atm, Spec_bound,
	     Sza, Zen, Azm, Pure_nadir, Rt->number_stokes(),
	     Use_first_order_scatt_calc,
	     Do_second_order, Number_stream,
	     Spectrum_spacing, Ps_mode, true),
  rt(Rt)
{
}

LRadRtBase::LRadRtBase
(const boost::shared_ptr<RtAtmosphere>& Atm,
 const SpectralBound& Spec_bound,
 const blitz::Array<double, 1>& Sza,
 const blitz::Array<double, 1>& Zen,
 const blitz::Array<double, 1>& Azm,
 bool Pure_nadir,
 int Number_stokes,
 bool Use_first_order_scatt_calc,
 bool Do_second_order,
 int Number_stream,
 double Spectrum_spacing,
 const LRadDriver::PsMode Ps_mode,
 bool Use_tms_correction)
  : atm_(Atm),
    ground_(Atm->ground()),
    use_first_order_scatt_calc(Use_first_order_scatt_calc),
    do_second_order(Do_second_order),
    sza(Sza.copy()), zen(Zen.copy()), azm(Azm.copy()),
    alt_spec_index_cache(-1)
{
  initialize(Spec_bound, Spectrum_spacing);

  // There seems to a bug in l_rad if nstokes = 4,
  // it will not compute the third stokes value,
  // so for now the maximum value is 3
  driver.reset(new LRadDriver(Number_stream,
			      std::min(Number_stokes, 3),
			      surface_type(),
			      Use_tms_correction,
			      Pure_nadir, Ps_mode));
}

//-----------------------------------------------------------------------
/// Constructor. Note that we can run in two modes, either we are
/// doing a polarization correction to an underlying multi-scatter RT
/// code, or we are just doing a single-scatter calculation in the
/// LRad alone. This constructor sets up for the single-scatter
/// calculation only, without a multi-scatter RT.
///
/// \param Stokes_coef The stokes coefficients to go from vector stokes
///    parameters to reflectance. This should be number_spectrometer() x 4.
/// \param Atm The RtAtmosphere to use.
/// \param Number_stream The number of streams to calculate with. Note
///    that this is the "half streams" more commonly used rather than
///    the "full streams" used in LIDORT 3.0.
/// \param Spec_bound Spectral window bounds for each spectrometer.
/// \param Sza Solar zenith angle. This is in degrees, and should be
///     in the range 0 to 90, and have size number_spectrometer()
/// \param Zen Zenith angle (degrees), in range 0 to 90, and have size
///      number_spectrometer()
/// \param Azm Azimuth angle (degrees), in range 0 to 360, and have size
///      number_spectrometer()
/// \param Pure_nadir If true then scene is purely nadir
/// \param Number_stokes Number of stokes coefficients to use.
/// \param Do_second_order If true, we do second order corrections
/// \param Spectrum_spacing The spectrum spacing
/// \param ps_mode Which pseudo spherical mode to use: 
///        REGULAR, ENHANCED, PLANE_PARALLEL, DETECT
//-----------------------------------------------------------------------

LRadRt::LRadRt(const boost::shared_ptr<StokesCoefficient>& Stokes_coef,
               const boost::shared_ptr<RtAtmosphere>& Atm,
               const SpectralBound& Spec_bound,
               const blitz::Array<double, 1>& Sza,
               const blitz::Array<double, 1>& Zen,
               const blitz::Array<double, 1>& Azm,
               bool Pure_nadir,
               int Number_stokes,
               bool Do_second_order,
               int Number_stream,
               double Spectrum_spacing,
               const LRadDriver::PsMode Ps_mode)
: RadiativeTransferSingleWn(Stokes_coef, Atm),
  LRadRtBase(Atm, Spec_bound, Sza, Zen, Azm, Pure_nadir, Number_stokes,
	     true, Do_second_order, Number_stream, Spectrum_spacing, Ps_mode,
	     false // Do not use TMS correction, no RT
	     )
{
}

void LRadRtBase::initialize
(const SpectralBound& Spec_bound, double Spectrum_spacing)
{
  if(Spec_bound.number_spectrometer() != atm_->number_spectrometer() ||
     sza.rows() != atm_->number_spectrometer() ||
     zen.rows() != atm_->number_spectrometer() ||
     azm.rows() != atm_->number_spectrometer())
    throw Exception("Spec_bound, Sza, Zen, and Atm all need to be size number_spectrometer()");

  for(int i = 0; i < atm_->number_spectrometer(); ++i) {
    range_check(sza(i), 0.0, 90.0);
    range_check(azm(i), 0.0, 360.0);
    range_check(zen(i), 0.0, 90.0);
  }

  if(!atm_->ground())
    throw Exception("LRadDriver cannot be used without a ground");

  for(int i = 0; i < atm_->number_spectrometer(); ++i) {
    // Some versions of the absco tables do not allow interpolation in
    // the wn direction. This means that only values with a wavenumber
    // a multiple of the grid spacing should be used.
    double t = Spec_bound.lower_bound(i, units::inv_cm).value;
    t = round(t / Spectrum_spacing) * Spectrum_spacing;
    wmin.push_back(t);
    
    t = Spec_bound.upper_bound(i, units::inv_cm).value;
    t = round(t / Spectrum_spacing) * Spectrum_spacing;
    wmax.push_back(t);
  }
  
  // Watch atmosphere for changes, so we clear cache if needed.
  atm_->add_observer(*this);
}

int LRadRtBase::surface_type() const
{
  if(!ground_)
    throw Exception("Need to have a ground to determine surface_type");
  SpurrBrdfType b = ground_->spurr_brdf_type();
  if(b == SpurrBrdfType::LAMBERTIAN)
    return BrdfType::LAMBERTIAN;
  if(b == SpurrBrdfType::COXMUNK)
    return BrdfType::COXMUNK;
  if(b == SpurrBrdfType::BREONVEG)
    return BrdfType::BPDFVEGN;
  if(b == SpurrBrdfType::BREONSOIL)
    return BrdfType::BPDFSOIL;
  Exception e;
  e << "Unrecognized SpurrBrdfType: " << b;
  throw e;
}

void LRadRtBase::setup_z_matrix_interpol
(const double wmin, const ArrayAd<double, 3>& pf_min, const double wmax,
 const ArrayAd<double, 3>& pf_max) const
{
  ArrayAd<double, 2> z_matrix_min(driver->z_matrix(pf_min));
  ArrayAd<double, 2> z_matrix_max(driver->z_matrix(pf_max));

  // Note that timings showed it was faster to have the value and
  // Jacobian interpolation done separately. This isn't dramatically
  // faster, but it is easy enough to do. The other approach would be
  // a single interpolation of an Array<AutoDerivative<double>, 2>.
  zmat_interpolate.reset(new
			 LinearInterpolate2Point<double, Array<double, 2> >
			 (wmin, z_matrix_min.value(),
			  wmax, z_matrix_max.value()));
  if(!z_matrix_min.is_constant())
    l_zmat_interpolate.reset(new
			     LinearInterpolate2Point<double, Array<double, 3> >
			     (wmin, z_matrix_min.jacobian(),
			      wmax, z_matrix_max.jacobian()));
}

//-----------------------------------------------------------------------
/// Update the altitude information. This can change the number of
/// layers if desired.
//-----------------------------------------------------------------------

void LRadRtBase::update_altitude(int spec_index) const
{
  if(spec_index == alt_spec_index_cache)
    return;

  alt_spec_index_cache = spec_index;

  Array<double, 1>
    alt(atm_->altitude(spec_index).convert(units::km).value.value());

  driver->setup_geometry(alt, sza(spec_index), zen(spec_index),
			 azm(spec_index));
  ArrayAd<double, 3> pf_min(atm_->phase_function_moments_wrt_rt(wmin[spec_index],
							       spec_index));
  ArrayAd<double, 3> pf_max(atm_->phase_function_moments_wrt_rt(wmax[spec_index],
							       spec_index));
   
  setup_z_matrix_interpol(wmin[spec_index], pf_min, wmax[spec_index], pf_max);
}

ArrayAd<double, 2> LRadRtBase::get_z_matrix
(const double Wn,
 int UNUSED(Spec_index),
 const boost::shared_ptr<OpticalProperties>& Opt_prop) const
{
  ArrayAd<double, 2> zmat;

  // If we aren't passed in intermediate variables to use by the LSI,
  // then we can use our interpolated value.
  if(!Opt_prop) {
    zmat.value().reference((*zmat_interpolate)(Wn));
    
    if(l_zmat_interpolate)
      zmat.jacobian().reference((*l_zmat_interpolate)(Wn));
  } else {
    // l_rad fortran code expects there to be 6 scattering moments, so
    // ensure that this number is passed through regardless of the source data
    ArrayAd<double, 3> pf(Opt_prop->total_phase_function_moments(-1, 6));
    zmat.reference(driver->z_matrix(pf));
  }

  return zmat;
}

void LRadRtBase::apply_jacobians
(double Wn, int Spec_index, ArrayAd<double, 1>& stokes,
 const Array<double, 3>& jac_atm, const Array<double, 2>& jac_surf,
 const boost::shared_ptr<OpticalProperties>& Opt_prop) const
{
  //-----------------------------------------------------------------------
  /// To speed up the calculation, the Atmosphere Jacobian was
  /// calculated relative to the RtAtmosphere "intermediate
  /// variables". The Surface Jacobian was calculated relative to the
  /// surface parameters. For both of these, calculate these relative to
  /// the state vector variables. Then sum the Atmosphere Jacobian over
  /// the layers and add in the surface Jacobian to give us the total
  /// Jacobian to the stokes with respect to the state vector.
  //-----------------------------------------------------------------------

  Array<double, 3> jac_iv(0, 0, 0);

  if(!Opt_prop) {
    Array<double, 3> inter_jac(atm_->intermediate_jacobian(Wn, Spec_index));

    if(inter_jac.cols() > 0)
      jac_iv.reference(inter_jac);
  } else if(Opt_prop->intermediate_jacobian().cols() > 0)
    jac_iv.reference(Opt_prop->intermediate_jacobian());
  if (stokes.number_variable() != jac_iv.depth())
    stokes.resize_number_variable(jac_iv.depth());

  ArrayAd<double, 1>
    surface_param(atm_->ground()->surface_parameter(Wn, Spec_index));

  for(int i = 0; i < stokes.rows(); ++i) {
    for(int j = 0; j < stokes.number_variable(); ++j) {
      double val = 0;
      if(jac_atm.depth() != 0)
	for(int m = 0; m < jac_iv.rows(); ++m)
	  for(int n = 0; n < jac_iv.cols(); ++n)
	    val += jac_atm(i, m, n) * jac_iv(m, n, j);
      if(!surface_param.is_constant())
	for(int m = 0; m < surface_param.jacobian().rows(); ++m)
	  val += jac_surf(i, m) * surface_param.jacobian()(m, j);
      stokes.jacobian()(i, j) = val;
    }
  }
}

//-----------------------------------------------------------------------
/// Calculate the stokes correction, which is pretty much all of
/// LRadRt expect for adding the actual underlying RT.
//-----------------------------------------------------------------------

blitz::Array<double, 1> LRadRtBase::stokes_corr
(double Wn, int Spec_index, bool have_rt,
 const boost::shared_ptr<OpticalProperties>& Opt_prop) const
{
  update_altitude(Spec_index);

  driver->setup_surface_params(atm_->ground()->surface_parameter(Wn, Spec_index).value());

  Array<double, 1> tau;
  Array<double, 1> omega;
  if(!Opt_prop) {
    tau.reference(atm_->optical_depth_wrt_rt(Wn, Spec_index).value());
    omega.reference(atm_->single_scattering_albedo_wrt_rt(Wn, Spec_index).value());
  } else {
    tau.reference(Opt_prop->total_optical_depth().value());
    omega.reference(Opt_prop->total_single_scattering_albedo().value());
  }
  Array<double, 3> pf(0, 0, 0);

  if(have_rt || do_second_order) {
    // -1 for numscat means get them all.
    int nscat = (do_second_order ? -1 : 1);
    int nmom = 2 * number_stream();

    if(!Opt_prop)
      pf.reference(atm_->phase_function_moments_wrt_rt(Wn, Spec_index, nmom, nscat).value());
    else
      pf.reference(Opt_prop->total_phase_function_moments(nmom, nscat).value());
  } 
  // For second order corrections, we need all the scattering
  // elements. However for first order we only need the first element,
  // and then only if we are calculating fscale. The simplest thing
  // would be to just get the full scattering matrix each time - but
  // this turns out to be a bit of a bottle neck so it is worth the
  // more complicated logic to only get what we need.
  Array<double, 2> zmat = get_z_matrix(Wn, Spec_index, Opt_prop).value();

  driver->setup_optical_inputs(tau, omega, pf, zmat);

  driver->clear_linear_inputs();
  
  // When not using first order calculation, zero out first stokes term,
  // This is disable in cases where first order of scattering comes from 
  // annother RT code or we are only concerned with the second order
  // correction
  if(use_first_order_scatt_calc) {
    driver->calculate_first_order();
  } else {
    driver->stokes()(0) = 0;
  }

  if(do_second_order)
    driver->calculate_second_order();

  // Copy value of stokes from driver so changes to it do not impact our answer
  // if we had used a reference
  Array<double, 1> stokes(driver->stokes().shape());
  stokes = driver->stokes();
  return stokes;
}

//-----------------------------------------------------------------------
/// Calculate the stokes correction, which is pretty much all of
/// LRadRt expect for adding the actual underlying RT.
//-----------------------------------------------------------------------

ArrayAd<double, 1> LRadRtBase::stokes_and_jacobian_corr
(double Wn, int Spec_index, bool have_rt,
 const boost::shared_ptr<OpticalProperties>& Opt_prop) const
{
  update_altitude(Spec_index);

  driver->setup_surface_params(atm_->ground()->surface_parameter(Wn, Spec_index).value());

  ArrayAd<double, 1> tau;
  ArrayAd<double, 1> omega;
  Array<double, 3> jac_iv(0, 0, 0);

  if(!Opt_prop) {
    tau.reference(atm_->optical_depth_wrt_rt(Wn, Spec_index));
    omega.reference(atm_->single_scattering_albedo_wrt_rt(Wn, Spec_index));
    Array<double, 3> inter_jac(atm_->intermediate_jacobian(Wn, Spec_index));

    if(inter_jac.cols() > 0)
      jac_iv.reference(inter_jac);
  } else {
    tau.reference(Opt_prop->total_optical_depth());
    omega.reference(Opt_prop->total_single_scattering_albedo());

    if(Opt_prop->intermediate_jacobian().cols() > 0)
      jac_iv.reference(Opt_prop->intermediate_jacobian());
  }

  ArrayAd<double, 3> pf(0, 0, 0, 0);

  if(have_rt || do_second_order) {
    // -1 for numscat means get them all.
    int nscat = (do_second_order ? -1 : 1);
    int nmom = 2 * number_stream();

    if(!Opt_prop)
      pf.reference(atm_->phase_function_moments_wrt_rt(Wn, Spec_index, nmom, nscat));
    else
      pf.reference(Opt_prop->total_phase_function_moments(nmom, nscat));
  } 

  // For second order corrections, we need all the scattering
  // elements. However for first order we only need the first element,
  // and then only if we are calculating fscale. The simplest thing
  // would be to just get the full scattering matrix each time - but
  // this turns out to be a bit of a bottle neck so it is worth the
  // more complicated logic to only get what we need.
  ArrayAd<double, 2> zmat = get_z_matrix(Wn, Spec_index, Opt_prop);

  driver->setup_optical_inputs(tau.value(), omega.value(), pf.value(), zmat.value());

  driver->setup_linear_inputs(tau, omega, pf, zmat);

  // When not using first order calculation, zero out first stokes term,
  // This is disable in cases where first order of scattering comes from 
  // annother RT code or we are only concerned with the second order
  // correction
  if(use_first_order_scatt_calc) {
    driver->calculate_first_order();
  } else {
    driver->stokes()(0) = 0;
    driver->surface_jacobian()(0, Range::all()) = 0;
    driver->atmospheric_jacobian()(0, Range::all(), Range::all()) = 0;
  }
  if(do_second_order)
    driver->calculate_second_order();
  
  // Calculate ArrayAd version of l_rad calculations that includes surface
  // and atmospheric jacobian components applied through the chain rule
  // using the intermediate jacobians
  ArrayAd<double, 1> stokes(driver->stokes().shape(), jac_iv.depth());
  stokes.value() = driver->stokes();
  apply_jacobians(Wn, Spec_index, stokes, driver->atmospheric_jacobian(),
                  driver->surface_jacobian(), Opt_prop);
  return stokes;
}

// See base class for description of this.
blitz::Array<double, 1> LRadRt::stokes_single_wn
(double Wn, int Spec_index,
 const boost::shared_ptr<OpticalProperties>& Opt_prop) const
{
  bool have_rt = (rt ? true : false);
  Array<double, 1> stokes = stokes_corr(Wn, Spec_index, have_rt, Opt_prop);
  if(rt) {
    Array<double, 1> t(rt->stokes_single_wn(Wn, Spec_index, Opt_prop));
    for(int i = 0; i < number_stokes(); ++i)
      stokes(i) = stokes(i) + t(i);
  }

  return stokes;
}

// See base class for description of this.
ArrayAd<double, 1> LRadRt::stokes_and_jacobian_single_wn
(double Wn, int Spec_index,
 const boost::shared_ptr<OpticalProperties>& Opt_prop) const
{
  bool have_rt = (rt ? true : false);
  ArrayAd<double, 1> stokes = stokes_and_jacobian_corr
    (Wn, Spec_index, have_rt, Opt_prop);
  
  // We either are correcting a multi-scatter RT code, or just doing a
  // single scatter alone.
  if(rt) {
    ArrayAd<double, 1> t(rt->stokes_and_jacobian_single_wn(Wn, Spec_index, Opt_prop));

    for(int i = 0; i < number_stokes(); ++i)
      stokes(i) = stokes(i) + t(i);
  }
  return stokes;
}

// See base class for description of this.
blitz::Array<double, 2> LRadRtFixedStokesCoefficient::stokes
(const SpectralDomain& Spec_domain, int Spec_index) const
{
  Array<double, 2> stokes = rt->stokes(Spec_domain, Spec_index);
  Array<double, 1> wn(Spec_domain.wavenumber());
  for(int i = 0; i < wn.rows(); ++i) {
    Array<double, 1> corr = stokes_corr(wn(i), Spec_index, true);
    for(int j = 0; j < number_stokes(); ++j)
      stokes(i,j) += corr(j);
  }
  return stokes;
}

// See base class for description of this.
ArrayAd<double, 2> LRadRtFixedStokesCoefficient::stokes_and_jacobian
(const SpectralDomain& Spec_domain, int Spec_index) const
{
  ArrayAd<double, 2> stokes = rt->stokes_and_jacobian(Spec_domain, Spec_index);
  Array<double, 1> wn(Spec_domain.wavenumber());
  for(int i = 0; i < wn.rows(); ++i) {
    ArrayAd<double, 1> corr = stokes_and_jacobian_corr(wn(i), Spec_index, true);
    for(int j = 0; j < number_stokes(); ++j)
      stokes(i,j) = stokes(i,j) + corr(j);
  }
  return stokes;
}

void LRadRtBase::print(std::ostream& Os, bool Short_form) const
{
    OstreamPad opad(Os, "  ");
    Os << "LRadRtBase:\n";

    opad << "Solar zenith angle:\n";
    opad << sza << "\n";
    opad.strict_sync();

    opad << "Zenith angle:\n";
    opad << zen << "\n";
    opad.strict_sync();

    opad << "Azimuth angle:\n";
    opad << azm << "\n";
    opad.strict_sync();

    opad << "Use first order scattering calculation: ";
    opad << (use_first_order_scatt_calc ? "True" : "False") << "\n\n";
    opad.strict_sync();

    opad << "Perform second order calculation: ";
    opad << (do_second_order ? "True" : "False") << "\n\n";
    opad.strict_sync();

    opad << "Z matrix interp wavenumber ranges:\n";
    for(int wn_idx = 0; wn_idx < (int) wmin.size(); wn_idx++) {
        opad << "  " << wmin[wn_idx] << ", " << wmax[wn_idx] << "\n";
    }
    opad << "\n";
    opad.strict_sync();

    // Output print from buffer class
    Os << "  Driver:\n";
    driver->print(opad, Short_form);
    opad << "\n";
    opad.strict_sync();
}

void LRadRt::print(std::ostream& Os, bool Short_form) const
{
    OstreamPad opad(Os, "  ");
    Os << "LRadRt:\n";
    LRadRtBase::print(opad, Short_form);
    opad.strict_sync();
    // Output print from parent class
    RadiativeTransferSingleWn::print(opad, Short_form);
    opad.strict_sync();
    if (rt) {
        Os << "  Radiative Transfer:\n";
        OstreamPad opad1(Os, "    ");
        rt->print(opad1, true);
        opad1 << "\n";
        opad1.strict_sync();
    }
}

void LRadRtFixedStokesCoefficient::print
(std::ostream& Os, bool Short_form) const
{
    OstreamPad opad(Os, "  ");
    Os << "LRadRt:\n";
    LRadRtBase::print(opad, Short_form);
    opad.strict_sync();
    // Output print from parent class
    RadiativeTransferFixedStokesCoefficient::print(opad, Short_form);
    opad.strict_sync();
    if (rt) {
        Os << "  Radiative Transfer:\n";
        OstreamPad opad1(Os, "    ");
        rt->print(opad1, true);
        opad1 << "\n";
        opad1.strict_sync();
    }
}
