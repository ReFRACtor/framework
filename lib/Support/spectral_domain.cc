#include "spectral_domain.h"
#include "fp_exception.h"
#include "fp_serialize_support.h"
#include "old_constant.h"
using namespace FullPhysics;
using namespace blitz;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void SpectralDomain::serialize(Archive& ar,
			      const unsigned int UNUSED(version))
{
  FP_GENERIC_BASE(SpectralDomain);
  ar & FP_NVP_(data) & FP_NVP_(sindex) & FP_NVP_(units);
}

FP_IMPLEMENT(SpectralDomain);
#endif

#ifdef HAVE_LUA
#include "register_lua.h"

std::string spec_domain_type_preference(const SpectralDomain& dom)
{
  if (dom.type_preference() == SpectralDomain::PREFER_WAVENUMBER)
    return "wavenumber";
  else if (dom.type_preference() == SpectralDomain::PREFER_WAVELENGTH)
    return "wavelength";
  else
    return "UNKNOWN";
}

REGISTER_LUA_CLASS(SpectralDomain)
.def("data", &SpectralDomain::data)
.def("units", &SpectralDomain::units)
.def("type_preference", &spec_domain_type_preference)
.def(luabind::constructor<const ArrayAd<double, 1>&, const Unit&>())
.def(luabind::constructor<const Array<double, 1>&>())
.def(luabind::constructor<const Array<double, 1>&, const Unit&>())
.def(luabind::constructor<const ArrayWithUnit<double, 1>&>())
.def(luabind::constructor<const ArrayAd<double, 1>&, const blitz::Array<int, 1>&,const Unit&>())
.def(luabind::constructor<const Array<double, 1>&, const blitz::Array<int, 1>&>())
.def(luabind::constructor<const Array<double, 1>&, const blitz::Array<int, 1>&, const Unit&>())
.def(luabind::constructor<const ArrayWithUnit<double, 1>&, const blitz::Array<int, 1>&>())
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Constructor. The units indicate if the data passed in is
/// wavenumber (1/Length) or wavelength (Length)
//-----------------------------------------------------------------------
SpectralDomain::SpectralDomain(const ArrayAd<double, 1>& Data,
                               const Unit& Units)
: data_(Data),
  units_(Units)
{
  if(!Units.is_commensurate(units::micron) &&
     !Units.is_commensurate(units::inv_cm)) {
    Exception e;
    e << "The units passed to SpectralDomain must be Length or Length^-1.\n"
      << "Units passed in: \n" << Units;
    throw e;
  }
}

//-----------------------------------------------------------------------
/// Constructor. The units indicate if the data passed in is
/// wavenumber (1/Length) or wavelength (Length)
//-----------------------------------------------------------------------
SpectralDomain::SpectralDomain(const ArrayAd<double, 1>& Data,
                               const blitz::Array<int, 1>& Sindex,
                               const Unit& Units)
: data_(Data),
  sindex_(Sindex),
  units_(Units)
{
  if(!Units.is_commensurate(units::micron) &&
     !Units.is_commensurate(units::inv_cm)) {
    Exception e;
    e << "The units passed to SpectralDomain must be Length or Length^-1.\n"
      << "Units passed in: \n" << Units;
    throw e;
  }
  if(data_.rows() != sindex_.rows())
    throw Exception("Data and Sindex must be the same size");
}

//-----------------------------------------------------------------------
/// Constructor. The units indicate if the data passed in is
/// wavenumber (1/Length) or wavelength (Length)
//-----------------------------------------------------------------------
 SpectralDomain::SpectralDomain(const blitz::Array<double, 1>& Data,
                               const Unit& Units)
: data_(Data),
  units_(Units)
{
  if(!Units.is_commensurate(units::micron) &&
     !Units.is_commensurate(units::inv_cm)) {
    Exception e;
    e << "The units passed to SpectralDomain must be Length or Length^-1.\n"
      << "Units passed in: \n" << Units;
    throw e;
  }
}

//-----------------------------------------------------------------------
/// Constructor. The units indicate if the data passed in is
/// wavenumber (1/Length) or wavelength (Length)
//-----------------------------------------------------------------------
 SpectralDomain::SpectralDomain(const blitz::Array<double, 1>& Data,
                                const blitz::Array<int, 1>& Sindex,
                                const Unit& Units)
: data_(Data),
  sindex_(Sindex),
  units_(Units)
{
  if(!Units.is_commensurate(units::micron) &&
     !Units.is_commensurate(units::inv_cm)) {
    Exception e;
    e << "The units passed to SpectralDomain must be Length or Length^-1.\n"
      << "Units passed in: \n" << Units;
    throw e;
  }
  if(data_.rows() != sindex_.rows())
    throw Exception("Data and Sindex must be the same size");
}

//-----------------------------------------------------------------------
/// Constructor. The units indicate if the data passed in is
/// wavenumber (1/Length) or wavelength (Length)
//-----------------------------------------------------------------------
SpectralDomain::SpectralDomain(const ArrayWithUnit<double, 1>& Data)
  : data_(Data.value),
    units_(Data.units)
{
  if(!units_.is_commensurate(units::micron) &&
     !units_.is_commensurate(units::inv_cm)) {
    Exception e;
    e << "The units passed to SpectralDomain must be Length or Length^-1.\n"
      << "Units passed in: \n" << units_;
    throw e;
  }
}

//-----------------------------------------------------------------------
/// Constructor. The units indicate if the data passed in is
/// wavenumber (1/Length) or wavelength (Length)
//-----------------------------------------------------------------------
SpectralDomain::SpectralDomain(const ArrayWithUnit<double, 1>& Data,
                               const blitz::Array<int, 1>& Sindex)
  : data_(Data.value),
    sindex_(Sindex),
    units_(Data.units)
{
  if(!units_.is_commensurate(units::micron) &&
     !units_.is_commensurate(units::inv_cm)) {
    Exception e;
    e << "The units passed to SpectralDomain must be Length or Length^-1.\n"
      << "Units passed in: \n" << units_;
    throw e;
  }
  if(data_.rows() != sindex_.rows())
    throw Exception("Data and Sindex must be the same size");
}

//-----------------------------------------------------------------------
/// Return data as the supplied the units.
//-----------------------------------------------------------------------

Array<double, 1> SpectralDomain::convert_wave(const Unit& Units) const
{
  if(units_.is_commensurate(Units))
    return Array<double, 1>(FullPhysics::conversion(units_, Units) * data_.value());
  else
    return Array<double, 1>(FullPhysics::conversion(1 / units_, Units) / data_.value());
}

//-----------------------------------------------------------------------
/// Return data as the supplied the units, including the Jacobian
//-----------------------------------------------------------------------

ArrayAd<double, 1> SpectralDomain::convert_wave_ad(const Unit& Units) const
{
  if(units_.is_commensurate(Units)) {
    blitz::Array<AutoDerivative<double>, 1> res(data_.rows());
    res = FullPhysics::conversion(units_, Units) * data_.to_array();
    return ArrayAd<double, 1>(res);
  } else {
    blitz::Array<AutoDerivative<double>, 1> res(data_.rows());
    res = FullPhysics::conversion(1 / units_, Units) / data_.to_array();
    return ArrayAd<double, 1>(res);
  }
}

//-----------------------------------------------------------------------
/// Return data as wavenumbers. You can optionally supply the units to
/// use.
/// Throws an error if the the optionally supplied units are not
/// commensurate with cm^-1
//-----------------------------------------------------------------------

Array<double, 1> SpectralDomain::wavenumber(const Unit& Units) const
{
  if (Units.is_commensurate(units::inv_cm))
    return convert_wave(Units);
  else {
    stringstream err_msg;
    err_msg << "Supplied units: " 
            << Units.name() 
            << " are not commensurate with target units: " 
            << units::inv_cm.name();
    throw Exception(err_msg.str());
  }
}

//-----------------------------------------------------------------------
/// Return data as wavelengths You can optionally supply the units to
/// use. 
/// Throws an error if the the optionally supplied units are not
/// commensurate with microns
//-----------------------------------------------------------------------

Array<double, 1> SpectralDomain::wavelength(const Unit& Units) const
{
  if (Units.is_commensurate(units::micron))
    return convert_wave(Units);
  else {
    stringstream err_msg;
    err_msg << "Supplied units: " 
            << Units.name() 
            << " are not commensurate with target units: " 
            << units::micron.name();
    throw Exception(err_msg.str());
  } 
}

//-----------------------------------------------------------------------
/// We may want to convert from photon number per second to radiance
/// units. This gives the factor to use in converting.
//-----------------------------------------------------------------------

ArrayWithUnit<double, 1> SpectralDomain::photon_to_radiance_factor() const
{
  using namespace units;
  // Desired output units. Note this is just a convenience, since we
  // typically give wavelength in microns and wavenumber in cm^-1. The
  // conversion could also be handled outside of this class.
  const Unit output_units = (W / inv_cm) / (ph / s / micron);
  const Unit wavenumber_unit = inv_cm;
  const DoubleWithUnit alpha = OldConstant::speed_of_light * OldConstant::planck;
  ArrayWithUnit<double, 1> res;
  res.value.reference
    (Array<double, 1>(alpha.value / wavenumber(wavenumber_unit)));
  res.units = alpha.units / wavenumber_unit / ph;
  return res.convert(output_units);
}

//-----------------------------------------------------------------------
/// Adds extra padding on either end of the grid to extend the grid
/// a certain amount of units. Tries to keep the same spacing of the
/// original grid. Can handle grids regardless if they are in increasing
/// or decreasing order.
/// 
/// The returned array has a sample_index value that maps the returned
/// SpectralDomain back to the original where padded indexes are marked
/// with -1.
//-----------------------------------------------------------------------

SpectralDomain SpectralDomain::add_padding(const DoubleWithUnit& padding)
{
    // Convert to nm but keep the same array ordering
    if(data_.rows() == 1) {
        throw Exception("SpectralDomain with one point can not be padded");
    }

    Array<double, 1> grid_in(data());

    // Continue spacing interval used at both ends
    double grid_spacing_left = abs(grid_in(1) - grid_in(0));
    double grid_spacing_right = abs(grid_in(grid_in.rows()-1) - grid_in(grid_in.rows()-2));

    double pad_amount = padding.convert(units_).value;

    int num_end_points = ceil(pad_amount / max(grid_spacing_left, grid_spacing_right));

    // Create the padded array and add the input grid in the middle
    ArrayAd<double, 1> grid_padded_ad(grid_in.rows() + 2*num_end_points, data_ad().number_variable());
    Array<double, 1> grid_padded(grid_padded_ad.value());

    Range orig_range = Range(num_end_points, grid_in.rows() + num_end_points - 1);

    grid_padded_ad.value()(orig_range) = data();
    if(!grid_padded_ad.is_constant()) {
        grid_padded_ad.jacobian() = 0;
        grid_padded_ad.jacobian()(orig_range, Range::all()) = data_ad().jacobian();
    }

    // Create indexing into original array for caller
    Array<int, 1> sample_indexes(grid_padded.rows());
    sample_indexes = -1;
    int orig_idx = 0;
    for(int samp_idx = num_end_points; samp_idx < grid_in.rows() + num_end_points; samp_idx++) {
        sample_indexes(samp_idx) = orig_idx++;
    }

    if(grid_in(1) > grid_in(0)) {
        // Grid in increasing order
        double curr_pad_val = grid_in(0) - num_end_points * grid_spacing_left;
        
        // Pad left of original grid
        for(int grid_idx = 0; grid_idx < num_end_points; grid_idx++) {
            grid_padded(grid_idx) = curr_pad_val;
            curr_pad_val += grid_spacing_left;
        }

        // Pad to the right of grid
        curr_pad_val = grid_in(grid_in.rows()-1) + grid_spacing_right;
        for(int grid_idx = grid_in.rows() + num_end_points; grid_idx < grid_in.rows() + num_end_points*2; grid_idx++) {
            grid_padded(grid_idx) = curr_pad_val;
            curr_pad_val += grid_spacing_right;
        }
    } else {
        // Grid in decreasing order 
        double curr_pad_val = grid_in(0) + num_end_points * grid_spacing_left;
        
        // Pad left of original grid
        for(int grid_idx = 0; grid_idx < num_end_points; grid_idx++) {
            grid_padded(grid_idx) = curr_pad_val;
            curr_pad_val -= grid_spacing_left;
        }

        // Pad to the right of grid
        curr_pad_val = grid_in(grid_in.rows()-1) - grid_spacing_right;
        for(int grid_idx = grid_in.rows() + num_end_points; grid_idx < grid_in.rows() + num_end_points*2; grid_idx++) {
            grid_padded(grid_idx) = curr_pad_val;
            curr_pad_val -= grid_spacing_right;
        }
    }

    return SpectralDomain(grid_padded_ad, sample_indexes, units_);
}

/// Clones object into a new copy

SpectralDomain SpectralDomain::clone() const
{ 
    if (sindex_.rows() > 0) {
        return SpectralDomain(data_.copy(), sindex_.copy(), units_); 
    } else {
        return SpectralDomain(data_.copy(), units_); 
    }
}
