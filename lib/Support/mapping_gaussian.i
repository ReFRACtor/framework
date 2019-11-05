%include "fp_common.i"

%{
#include "mapping_gaussian.h"
%}

%base_import(mapping)
%import "array_ad.i"
%import "pressure.i"
%fp_shared_ptr(FullPhysics::MappingGaussian);


namespace FullPhysics {

%feature("notabstract") MappingGaussian;

class MappingGaussian : public Mapping  {
public:
    MappingGaussian(const boost::shared_ptr<Pressure>& in_press, bool Linear_AOD);
    virtual ~MappingGaussian();
    virtual const ArrayAd<double, 1> fm_view(ArrayAd<double, 1> const& updated_coeff) const;
    virtual const ArrayAd<double, 1> retrieval_init(ArrayAd<double, 1> const& initial_coeff) const;
    virtual std::string name() const;
    virtual boost::shared_ptr<Mapping> clone() const;
    virtual bool is_linear_total() const;
    virtual AutoDerivative<double> total_optical_depth(ArrayAd<double, 1> component) const;
};
const double MappingGaussian::min_desired = 1e-9;
}
