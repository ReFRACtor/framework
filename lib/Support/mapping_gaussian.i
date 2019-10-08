%include "fp_common.i"

%{
#include "mapping_gaussian.h"
%}

%base_import(mapping)
%import "array_ad.i"
%import "pressure.i"
%fp_shared_ptr(FullPhysics::MappingGaussian);


namespace FullPhysics {
class MappingGaussian : public Mapping  {
public:
    MappingGaussian(const boost::shared_ptr<Pressure>& in_press, bool Linear_AOD);
    const ArrayAd<double, 1> fm_view(ArrayAd<double, 1> const& updated_coeff) const;
    const ArrayAd<double, 1> retrieval_init(ArrayAd<double, 1> const& initial_coeff) const;
    std::string name() const;
    boost::shared_ptr<Mapping> clone() const;
    bool is_linear_total() const;
    AutoDerivative<double> total_optical_depth(ArrayAd<double, 1> component) const;
    ~MappingGaussian();
};
const double MappingGaussian::min_desired = 1e-9;
}
