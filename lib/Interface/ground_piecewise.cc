#include "ground_piecewise.h"
#include "fp_serialize_support.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void GroundPiecewise::serialize(Archive & ar,
			const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(GroundImpBase)
    & FP_NVP_(spectral_points);
  // Skip saving this, instead we recreate it.
  //    & FP_NVP(ground_interp);
}

template<class Archive>
void GroundPiecewise::save(Archive & UNUSED(ar),
			  const unsigned int UNUSED(version)) const
{
  // Nothing more to do
}
template<class Archive>
void GroundPiecewise::load(Archive & UNUSED(ar),
			  const unsigned int UNUSED(version))
{
  // This updates ground_interp.
  update_sub_state_hook();
}


FP_IMPLEMENT(GroundPiecewise);
#endif

GroundPiecewise::GroundPiecewise(const ArrayWithUnit<double, 1>& spectral_points,
                                 const blitz::Array<double, 1>& point_values,
                                 const boost::shared_ptr<StateMapping>& mapping)
{
  spectral_points_ = spectral_points;
  init(point_values, mapping);
  update_sub_state_hook();
}


/// Compute surface parameter array used by radiative transfer
ArrayAd<double, 1> GroundPiecewise::surface_parameter
(const double wn, const int UNUSED(spec_index)) const
{
  AutoDerivative<double> wn_value = value_at_point(DoubleWithUnit(wn, units::inv_cm));
  ArrayAd<double, 1> spars(1, wn_value.number_variable());
  spars(0) = wn_value;
  return spars;
}

/// Return value by querying with an spectral point in arbitrary units
const AutoDerivative<double> GroundPiecewise::value_at_point(const DoubleWithUnit wave_point) const
{
  double point_val = wave_point.convert_wave(spectral_points_.units).value;
  return (*ground_interp)(point_val);
}

void GroundPiecewise::update_sub_state_hook()
{
  ArrayAd<double, 1> mapped_state_coeff = mapping->mapped_state(coeff);

  std::vector<AutoDerivative<double> > points_list;
  std::vector<AutoDerivative<double> > coeff_list;
  for(int idx = 0; idx < coefficient().rows(); idx++) {
    points_list.push_back(spectral_points_.value(idx));
    coeff_list.push_back(mapped_state_coeff(idx));
  }

  typedef LinearInterpolate<AutoDerivative<double>, AutoDerivative<double> > interp_type;
  ground_interp.reset( new interp_type(points_list.begin(), points_list.end(), coeff_list.begin()) );
}
