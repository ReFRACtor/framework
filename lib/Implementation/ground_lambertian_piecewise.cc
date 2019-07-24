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
    return boost::shared_ptr<Ground>(new GroundLambertianPiecewise(ArrayWithUnit<double, 1>(wavenumbers, units::inv_cm), coefficient().value(), used_flag));
}

std::string GroundLambertianPiecewise::sub_state_identifier() const {
    return "ground/lambertian_piecewise";
}

std::string GroundLambertianPiecewise::state_vector_name_i(int i) const
{
    return "Ground Lambertian Piecewise at wavenumber " + boost::lexical_cast<std::string>(wavenumbers(i));
}

void GroundLambertianPiecewise::print(std::ostream& Os) const
{
    OstreamPad opad(Os, "    ");
    Os << "GroundLambertianPiecewise:" << std::endl;

    opad << wavenumbers.rows() << " lambertian coefficients (" << wavenumbers(0) << "-" << wavenumbers(wavenumbers.rows()-1) << " cm^-1)" << std::endl;

    opad.strict_sync();
}

std::string GroundLambertianPiecewise::desc() const {
    return "GroundLambertianPiecewise";
}
