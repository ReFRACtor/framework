#include "ground_lambertian_piecewise.h"
#include "ostream_pad.h"
#include <boost/lexical_cast.hpp>

using namespace FullPhysics;
using namespace blitz;

GroundLambertianPiecewise::GroundLambertianPiecewise(const ArrayWithUnit<double, 1>& spectral_points,
                                                     const blitz::Array<double, 1>& lambertian_values,
                                                     const blitz::Array<bool, 1>& retrieval_flag)
: GroundPiecewise(spectral_points, lambertian_values, retrieval_flag)
{
}

boost::shared_ptr<Ground> GroundLambertianPiecewise::clone() const
{
    return boost::shared_ptr<Ground>(new GroundLambertianPiecewise(spectral_points_, coefficient().value(), used_flag));
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

std::string GroundLambertianPiecewise::desc() const {
    return "GroundLambertianPiecewise";
}
