#include "atmosphere_standard.h"
#include "fp_serialize_support.h"

#include <boost/foreach.hpp>
#include <cmath>

#include "constant.h"
#include "linear_algebra.h"
#include "ground.h"
#include "planck.h"
#include "altitude.h"
#include "aerosol_optical.h"
#include "ostream_pad.h"
#include "logger.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void AtmosphereStandard::serialize(Archive & ar,
				   const unsigned int version)
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(RtAtmosphere)
    & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ObserverAerosol)
    & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ObserverPressure);
  // Older version didn't ObserverAbsorber. 
  if(version > 0)
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ObserverAbsorber);
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(CacheInvalidatedObserver)
    & FP_NVP(absorber)
    & FP_NVP(pressure)
    & FP_NVP(temperature)
    & FP_NVP(aerosol)
    & FP_NVP(rh)
    & FP_NVP(ground_ptr)
    & FP_NVP(rayleigh)
    & FP_NVP(surface_temp)
    & FP_NVP(constant)
    & FP_NVP(alt)
    & FP_NVP(sv_jac_size)
    & FP_NVP(nlay)
    & FP_NVP(can_cache_channel);
}

FP_IMPLEMENT(AtmosphereStandard);
#endif

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(AtmosphereStandard, RtAtmosphere)
.def(luabind::constructor<const boost::shared_ptr<Absorber>&,
     const boost::shared_ptr<Pressure>&,
     const boost::shared_ptr<Temperature>&,
     const boost::shared_ptr<Rayleigh>&,
     const boost::shared_ptr<Aerosol>&,
     const boost::shared_ptr<RelativeHumidity>&,
     const boost::shared_ptr<Ground>&,
     const std::vector<boost::shared_ptr<Altitude> >&,
     const boost::shared_ptr<Constant>&>())
.def(luabind::constructor<const boost::shared_ptr<Absorber>&,
     const boost::shared_ptr<Pressure>&,
     const boost::shared_ptr<Temperature>&,
     const boost::shared_ptr<Rayleigh>&,
     const boost::shared_ptr<RelativeHumidity>&,
     const boost::shared_ptr<Ground>&,
     const std::vector<boost::shared_ptr<Altitude> >&,
     const boost::shared_ptr<Constant>&>())
.def(luabind::constructor<const boost::shared_ptr<Absorber>&,
     const boost::shared_ptr<Pressure>&,
     const boost::shared_ptr<Temperature>&,
     const boost::shared_ptr<Rayleigh>&,
     const boost::shared_ptr<Aerosol>&,
     const boost::shared_ptr<RelativeHumidity>&,
     const std::vector<boost::shared_ptr<Altitude> >&,
     const boost::shared_ptr<Constant>&>())
.def(luabind::constructor<const boost::shared_ptr<Absorber>&,
     const boost::shared_ptr<Pressure>&,
     const boost::shared_ptr<Temperature>&,
     const boost::shared_ptr<Rayleigh>&,
     const boost::shared_ptr<RelativeHumidity>&,
     const std::vector<boost::shared_ptr<Altitude> >&,
     const boost::shared_ptr<Constant>&>())
.def("pressure", &AtmosphereStandard::pressure_ptr)
.def("temperature", &AtmosphereStandard::temperature_ptr)
.def("ground", &AtmosphereStandard::ground)
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Create an an Atmosphere with all available components:
///
/// Required:
/// * Pressure
/// * Temperature
/// * Rayleigh
/// * Relative Humdity
/// * Altitude
///
/// Optional:
/// * Aerosol
/// * Ground
/// * Surface temperature
//-----------------------------------------------------------------------

AtmosphereStandard::AtmosphereStandard(const boost::shared_ptr<Absorber>& absorberv,
                                       const boost::shared_ptr<Pressure>& pressurev,
                                       const boost::shared_ptr<Temperature>& temperaturev,
                                       const boost::shared_ptr<Rayleigh>& rayleighv,
                                       const boost::shared_ptr<Aerosol>& aerosolv,
                                       const boost::shared_ptr<RelativeHumidity>& rhv,
                                       const boost::shared_ptr<Ground>& groundv,
                                       const boost::shared_ptr<SurfaceTemperature>& surface_tempv,
                                       const std::vector<boost::shared_ptr<Altitude> >& altv,
                                       const boost::shared_ptr<Constant>& C)
: absorber(absorberv), pressure(pressurev), temperature(temperaturev), 
  aerosol(aerosolv), rh(rhv), ground_ptr(groundv), rayleigh(rayleighv),
  surface_temp(surface_tempv),
  constant(C), alt(altv),
  sv_jac_size(0),
  nlay(-1)
{
    initialize();
}

//-----------------------------------------------------------------------
/// Create an Atmosphere with required components and all optional
/// components except for surface temperature.
//-----------------------------------------------------------------------

AtmosphereStandard::AtmosphereStandard(const boost::shared_ptr<Absorber>& absorberv,
                                       const boost::shared_ptr<Pressure>& pressurev,
                                       const boost::shared_ptr<Temperature>& temperaturev,
                                       const boost::shared_ptr<Rayleigh>& rayleighv,
                                       const boost::shared_ptr<Aerosol>& aerosolv,
                                       const boost::shared_ptr<RelativeHumidity>& rhv,
                                       const boost::shared_ptr<Ground>& groundv,
                                       const std::vector<boost::shared_ptr<Altitude> >& altv,
                                       const boost::shared_ptr<Constant>& C)
: absorber(absorberv), pressure(pressurev), temperature(temperaturev), 
  aerosol(aerosolv), rh(rhv), ground_ptr(groundv), rayleigh(rayleighv),
  constant(C), alt(altv),
  sv_jac_size(0),
  nlay(-1)
{
    initialize();
}

//-----------------------------------------------------------------------
/// Create an Atmosphere with required components and all optional
/// components except for ground and surface temperature.
///-----------------------------------------------------------------------

AtmosphereStandard::AtmosphereStandard(const boost::shared_ptr<Absorber>& absorberv,
                                       const boost::shared_ptr<Pressure>& pressurev,
                                       const boost::shared_ptr<Temperature>& temperaturev,
                                       const boost::shared_ptr<Rayleigh>& rayleighv,
                                       const boost::shared_ptr<Aerosol>& aerosolv,
                                       const boost::shared_ptr<RelativeHumidity>& rhv,
                                       const std::vector<boost::shared_ptr<Altitude> >& altv,
                                       const boost::shared_ptr<Constant>& C)
: absorber(absorberv), pressure(pressurev), temperature(temperaturev), 
  aerosol(aerosolv), rh(rhv), rayleigh(rayleighv),
  constant(C), alt(altv),
  sv_jac_size(0),
  nlay(-1)
{
    initialize();
}

//-----------------------------------------------------------------------
/// Create an Atmosphere with required components and all optional
/// components except for aerosol
///-----------------------------------------------------------------------

AtmosphereStandard::AtmosphereStandard(const boost::shared_ptr<Absorber>& absorberv,
                                       const boost::shared_ptr<Pressure>& pressurev,
                                       const boost::shared_ptr<Temperature>& temperaturev,
                                       const boost::shared_ptr<Rayleigh>& rayleighv,
                                       const boost::shared_ptr<RelativeHumidity>& rhv,
                                       const boost::shared_ptr<Ground>& groundv,
                                       const boost::shared_ptr<SurfaceTemperature>& surface_tempv,
                                       const std::vector<boost::shared_ptr<Altitude> >& altv,
                                       const boost::shared_ptr<Constant>& C)
: absorber(absorberv), pressure(pressurev), temperature(temperaturev), 
  rh(rhv), ground_ptr(groundv), rayleigh(rayleighv),
  surface_temp(surface_tempv),
  constant(C), alt(altv),
  sv_jac_size(0),
  nlay(-1)
{
    initialize();
}

//-----------------------------------------------------------------------
/// Create an Atmosphere with required components and all optional
/// components except for aerosol and surface temperature
///-----------------------------------------------------------------------

AtmosphereStandard::AtmosphereStandard(const boost::shared_ptr<Absorber>& absorberv,
                                       const boost::shared_ptr<Pressure>& pressurev,
                                       const boost::shared_ptr<Temperature>& temperaturev,
                                       const boost::shared_ptr<Rayleigh>& rayleighv,
                                       const boost::shared_ptr<RelativeHumidity>& rhv,
                                       const boost::shared_ptr<Ground>& groundv,
                                       const std::vector<boost::shared_ptr<Altitude> >& altv,
                                       const boost::shared_ptr<Constant>& C)
: absorber(absorberv), pressure(pressurev), temperature(temperaturev), 
  rh(rhv), ground_ptr(groundv),rayleigh(rayleighv),
  constant(C), alt(altv),
  sv_jac_size(0),
  nlay(-1)
{
    initialize();
}

//-----------------------------------------------------------------------
/// Create an Atmosphere with required components and all optional
/// components except for ground, aerosol and surface temperature
//-----------------------------------------------------------------------

AtmosphereStandard::AtmosphereStandard(const boost::shared_ptr<Absorber>& absorberv,
                                       const boost::shared_ptr<Pressure>& pressurev,
                                       const boost::shared_ptr<Temperature>& temperaturev,
                                       const boost::shared_ptr<Rayleigh>& rayleighv,
                                       const boost::shared_ptr<RelativeHumidity>& rhv,
                                       const std::vector<boost::shared_ptr<Altitude> >& altv,
                                       const boost::shared_ptr<Constant>& C)
: absorber(absorberv), pressure(pressurev), temperature(temperaturev), 
  rh(rhv), rayleigh(rayleighv), constant(C), alt(altv), sv_jac_size(0), nlay(-1)
{
    initialize();
}

void AtmosphereStandard::initialize()
{
  if(!absorber)
    throw Exception("Absorber is not allowed to be null in AtmosphereStandard");
  absorber->add_observer(*this);
  absorber->add_cache_invalidated_observer(*this);
  if(!pressure)
    throw Exception("Pressure is not allowed to be null in AtmosphereStandard");
  pressure->add_cache_invalidated_observer(*this);
  if(!temperature)
    throw Exception("Temperature is not allowed to be null in AtmosphereStandard");
  temperature->add_cache_invalidated_observer(*this);
    
  BOOST_FOREACH(const boost::shared_ptr<Altitude>& a, alt) {
      if(!a)
        throw Exception("Altitude is not allowed to be null in AtmosphereStandard");
      a->add_cache_invalidated_observer(*this);
  }

  if(aerosol) {
    aerosol->add_observer(*this);
    aerosol->add_cache_invalidated_observer(*this);
  }
  pressure->add_observer(*this);

  // Use a map caches, to allow the reuse of values needed outside of the
  // main spectral grid loops
  column_od_cache.reset( new ArrayAdMapCache<double, double, 1>() );
  total_od_cache.reset( new ArrayAdMapCache<double, double, 1>() );

  // Initially can not use the above two caches unless a StateVector gets attached
  // to handle invalidating these caches
  can_cache_channel = false;

  // Make sure caches are in a consistent state, zero out so that unintialized values
  // do not cause the debugging stats to print prematurely
  n_cache_calls = 0;
  invalidate_local_cache();
}

AtmosphereStandard::AtmosphereStandard()
{
  // Set up cache as empty, but don't do any other initialization.
 
  // Use a map caches, to allow the reuse of values needed outside of the
  // main spectral grid loops
  column_od_cache.reset( new ArrayAdMapCache<double, double, 1>() );
  total_od_cache.reset( new ArrayAdMapCache<double, double, 1>() );

  // Initially can not use the above two caches unless a
  /// StateVector gets attached
  // to handle invalidating these caches
  can_cache_channel = false;

  // Make sure caches are in a consistent state
  n_cache_calls = 0;
  invalidate_local_cache();
}

//-----------------------------------------------------------------------
/// Changes the aerosol class used. Notifies the relevant obsevers
/// and invalidates the cache
//-----------------------------------------------------------------------

void AtmosphereStandard::set_aerosol(boost::shared_ptr<Aerosol>& new_aerosol, StateVector& Sv)
{
    // Remove observers from old aerosol instance
    Sv.remove_observer(*aerosol);

    // Switch to new instance
    aerosol = new_aerosol;

    // Register observers
    aerosol->add_observer(*this);

    // Run notify update so gradient is set up correctly
    aerosol->notify_update(Sv);

    // Invalidate caches
    invalidate_local_cache();
}

void AtmosphereStandard::reset_timer()
{
    RtAtmosphere::reset_timer();
    Aerosol::timer.reset_elapsed();
    Absorber::timer.reset_elapsed();
}

std::string AtmosphereStandard::timer_info() const
{
    std::ostringstream os;
    os << RtAtmosphere::timer_info() << "\n"
       << "   " << Absorber::timer << "\n"
       << "   " << Aerosol::timer << "\n";
    return os.str();
}

//-----------------------------------------------------------------------
/// For performance, we cache some data as we calculate it. This
/// becomes stale when the aerosol is changed, so we observe aerosol
/// and mark the cache when it changes.
//-----------------------------------------------------------------------

void AtmosphereStandard::notify_update(const Aerosol& UNUSED(A))
{
    invalidate_local_cache();
    notify_update_do(*this);
}

//-----------------------------------------------------------------------
/// This handles filling the various cached variables. We first check
/// to see if these values are already cached, if so then we skip the
/// calculation.
/// Returns true if cache is filled, otherwise returns false
//-----------------------------------------------------------------------

bool AtmosphereStandard::fill_cache(double wn, int spec_index) const
{
    // If spectrometer changes then reset caches 
    if(spec_index != spec_index_tau_cache) {
        invalidate_local_cache();
        spec_index_tau_cache = spec_index;
    }

    n_cache_calls += 1;

    typename std::map<double, boost::shared_ptr<OpticalPropertiesWrtRt> >::const_iterator value_it;
    value_it = opt_prop_cache.find(wn);
    if (value_it != opt_prop_cache.end()) {

        opt_prop = value_it->second;

        n_cache_hits += 1;

        return true;

    } else {

        FunctionTimer ft(timer.function_timer());

        opt_prop.reset(new OpticalPropertiesWrtRt());
        opt_prop->initialize(DoubleWithUnit(wn, units::inv_cm), spec_index, absorber, rayleigh, aerosol, sv_jac_size);

        opt_prop_cache.insert(std::pair<double, boost::shared_ptr<OpticalPropertiesWrtRt> >(wn, opt_prop));

        n_cache_miss += 1;

    }

    return false;
}

//-----------------------------------------------------------------------
/// Reset cache flags so that cached values are recomputed upon next
/// query
//-----------------------------------------------------------------------

void AtmosphereStandard::invalidate_local_cache() const
{
    if (n_cache_calls > 0) {
        Logger::info() << "Resetting optical properties cache: "
            << n_cache_calls << " calls, " << n_cache_hits << " cache hits, " 
            << n_cache_miss << " cache misses.\n";
    }

    spec_index_tau_cache = -1;

    opt_prop_cache.clear();
    n_cache_calls = 0;
    n_cache_hits = 0;
    n_cache_miss = 0;

    column_od_cache->clear();
    total_od_cache->clear();
}

//
// See bass class for description
ArrayAdWithUnit<double, 1> AtmosphereStandard::altitude(int spec_index) const
{
    range_check(spec_index, 0, number_spectrometer());
    ArrayAdWithUnit<double, 1> p = pressure->pressure_grid();
    Array<AutoDerivative<double>, 1> res(p.rows());
    Unit u = alt[spec_index]->altitude(p(0)).units;

    for(int i = 0; i < res.rows(); ++i) {
        res(i) = alt[spec_index]->altitude(p(i)).convert(u).value;
    }

    return ArrayAdWithUnit<double, 1>(ArrayAd<double, 1>(res), u);
}

ArrayAd<double, 1> AtmosphereStandard::optical_depth_wrt_rt(double wn, int spec_index) const
{
    if (!can_cache_channel || !total_od_cache->is_valid(wn)) {
        fill_cache(wn, spec_index);

        total_od_cache->insert(wn, opt_prop->total_optical_depth());
    }

    return (*total_od_cache)[wn];
}

AutoDerivative<double> AtmosphereStandard::column_optical_depth(double wn, int spec_index, const std::string& Gas_name) const
{

    if (!can_cache_channel || !column_od_cache->is_valid(wn)) {
        firstIndex i1; secondIndex i2; thirdIndex i3; fourthIndex i4;

        // It is easier to go back to the absorber object to compute the jacobians correctly since
        // the jacobians in opt_prop are wrt to total gas optical depth
        ArrayAd<double, 2> taug_i = absorber->optical_depth_each_layer(wn, spec_index);
        ArrayAd<double, 1> totaltaug(taug_i.cols(), taug_i.number_variable());

        totaltaug.value() = sum(taug_i.value()(i2, i1), i2);

        if(!taug_i.is_constant()) {
            totaltaug.jacobian() = sum(taug_i.jacobian()(i3, i1, i2), i3);
        } else {
            totaltaug.jacobian() = 0;
        }

        column_od_cache->insert(wn, totaltaug);
    }

    return (*column_od_cache)[wn](absorber->gas_index(Gas_name));
}

//-----------------------------------------------------------------------
/// The atmospheric thermal blackbody values per level.
//-----------------------------------------------------------------------

ArrayAd<double, 1> AtmosphereStandard::atmosphere_blackbody(double wn, int UNUSED(spec_index)) const
{
    ArrayAdWithUnit<double, 1> temp_grid = temperature->temperature_grid(*pressure);

    ArrayAd<double, 1> black_body(temp_grid.rows(), temp_grid.number_variable());

    for(int lev_idx = 0; lev_idx < temp_grid.rows(); lev_idx++) {
        double temp_K = temp_grid(lev_idx).convert(Unit("K")).value.value();
        black_body(lev_idx) = planck(wn, temp_K, temp_grid(lev_idx).value.gradient());
    }

    return black_body;
}

//-----------------------------------------------------------------------
/// The surface thermal blackbody.
//-----------------------------------------------------------------------

AutoDerivative<double> AtmosphereStandard::surface_blackbody(double wn, int spec_index) const
{
    if (!surface_temp) {
        Exception error;
        error << "Cannot compute surface_blackbody, atmosphere was constructed without a surface tempererature object.";
        throw error;
    }

    double surf_temp_K = surface_temp->surface_temperature(spec_index).convert(Unit("K")).value.value();
    Array<double, 1> surf_temp_grad = surface_temp->surface_temperature(spec_index).value.gradient();
    return planck(wn, surf_temp_K, surf_temp_grad);
}

//-----------------------------------------------------------------------
/// This clones a Atmosphere object. This is a deep copy, all of the
/// objects that are part of this are cloned also (e.g., Pressure,
/// Temperature).
///
/// This cloned copy will *not* be attached to any StateVector, nor
/// will any Observer<Atmosphere> objects be attached (although the
/// original object is unchanged). You can attach the clone to any
/// objects you wish to.
///
/// This property is particularly useful to "freeze" the state. For
/// example, if the StateVector was set the apriori state and the
/// Atmosphere attached to the StateVector is cloned, then the cloned
/// version will continue to be the "Apriori atmosphere", even if the
/// StateVector is subsequently changed thus updating the original
/// object.
//-----------------------------------------------------------------------

boost::shared_ptr<AtmosphereStandard> AtmosphereStandard::clone() const
{
    boost::shared_ptr<Pressure> pressure_clone = pressure->clone();
    boost::shared_ptr<Temperature> temperature_clone =
        temperature->clone();
    boost::shared_ptr<Rayleigh> rayleigh_clone = rayleigh->clone();

    boost::shared_ptr<Ground> ground_clone;

    if(ground_ptr) {
        ground_clone = ground_ptr->clone();
    }

    std::vector<boost::shared_ptr<Altitude> > alt_clone;
    BOOST_FOREACH(const boost::shared_ptr<Altitude>& a, alt)
    alt_clone.push_back(a->clone());
    boost::shared_ptr<Absorber> absorber_clone =
        absorber->clone();
    boost::shared_ptr<RelativeHumidity> rh_clone =
        rh->clone();
    boost::shared_ptr<Aerosol> aerosol_clone;

    if(aerosol) {
        aerosol_clone = aerosol->clone();
    }

    boost::shared_ptr<AtmosphereStandard> res
    (new AtmosphereStandard(absorber_clone, pressure_clone, temperature_clone, rayleigh_clone,
                            aerosol_clone, rh_clone, ground_clone, alt_clone,
                            constant));
    return res;
}

void AtmosphereStandard::print(std::ostream& Os) const
{
    Os << "AtmosphereStandard:\n";
    OstreamPad opad(Os, "    ");
    Os << "  Constant:\n";
    opad << *constant << "\n";
    opad.strict_sync();
    Os << "  Absorber:\n";
    opad << *absorber << "\n";
    opad.strict_sync();
    Os << "  Pressure:\n";
    opad << *pressure << "\n";
    opad.strict_sync();
    Os << "  Temperature:\n";
    opad << *temperature << "\n";
    opad.strict_sync();
    Os << "  Aerosol:\n";

    if(aerosol) {
        opad << *aerosol << "\n";
    }
    else {
        opad << "Rayleigh only, no aerosols\n";
    }

    opad.strict_sync();
    Os << "  Ground:\n";

    if(ground_ptr) {
        opad << *ground_ptr << "\n";
    }
    else {
        opad << "No ground\n";
    }

    opad.strict_sync();

    for(int i = 0; i < (int) alt.size(); ++i) {
        Os << "  Altitude[" << i << "]:\n";
        opad << *(alt[i]) << "\n";
        opad.strict_sync();
    }
}

//-----------------------------------------------------------------------
/// For unit test purposes, it is useful to be able to clone or create
/// a new atmosphere class then attach all its nested children to the
/// state vector in the historic order and in the manner done by the
/// Lua configuration
//-----------------------------------------------------------------------

void AtmosphereStandard::attach_children_to_sv(StateVector& statev)
{
    statev.add_observer(*this);

    // Need to reattach the cloned atmosphere components to the SV one by one as done in the Lua config and in the same order
    statev.add_observer(*this->absorber_ptr());

    for (int spec_idx = 0; spec_idx < this->absorber_ptr()->number_species(); spec_idx++) {
        std::string gas_name = this->absorber_ptr()->gas_name(spec_idx);
        statev.add_observer(*this->absorber_ptr()->absorber_vmr(gas_name));
    }

    statev.add_observer(*this->pressure_ptr());
    statev.add_observer(*this->temperature_ptr());

    if(this->aerosol_ptr()) {
        statev.add_observer(*this->aerosol_ptr());

        for (int aer_idx = 0; aer_idx < this->aerosol_ptr()->number_particle(); aer_idx++) {
            boost::shared_ptr<AerosolOptical> aer_optical(boost::dynamic_pointer_cast<AerosolOptical>(this->aerosol_ptr()));

            if (aer_optical) {
                statev.add_observer(*aer_optical->aerosol_extinction(aer_idx));
                statev.add_observer(*aer_optical->aerosol_property(aer_idx));
            }
        }
    }

    if(this->ground()) {
        statev.add_observer(*this->ground());
    }
}
