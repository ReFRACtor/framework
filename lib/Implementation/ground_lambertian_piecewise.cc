#include "ground_lambertian_piecewise.h"
#include "fp_serialize_support.h"
#include "ostream_pad.h"
#include <boost/lexical_cast.hpp>

using namespace FullPhysics;
using namespace blitz;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void GroundLambertianPiecewise::serialize(Archive & ar,
			const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(GroundPiecewise);
}

FP_IMPLEMENT(GroundLambertianPiecewise);
#endif

GroundLambertianPiecewise::GroundLambertianPiecewise(const ArrayWithUnit<double, 1>& spectral_points,
                                                     const blitz::Array<double, 1>& lambertian_values,
                                                     const boost::shared_ptr<StateMapping>& mapping)
: GroundPiecewise(spectral_points, lambertian_values, mapping)
{
}

boost::shared_ptr<Ground> GroundLambertianPiecewise::clone() const
{
    return boost::shared_ptr<Ground>(new GroundLambertianPiecewise(spectral_points_, coefficient().value()));
}

std::string GroundLambertianPiecewise::sub_state_identifier() const {
    return "ground/lambertian_piecewise";
}

std::string GroundLambertianPiecewise::state_vector_name_i(int i) const
{
    return "Ground Lambertian Piecewise at " + boost::lexical_cast<std::string>(spectral_points_.value(i)) + " " + spectral_points_.units.name();
}

void GroundLambertianPiecewise::print(std::ostream& Os) const
{
    OstreamPad opad(Os, "    ");
    Os << "GroundLambertianPiecewise:" << std::endl;

    Array<double, 1> point_vals = spectral_points_.value;
    opad << point_vals.rows() << " lambertian coefficients (" << point_vals(0) << "-" << point_vals(point_vals.rows()-1) << " " << spectral_points_.units.name() << ")" << std::endl;

    opad.strict_sync();
}
