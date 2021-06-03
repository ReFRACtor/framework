%include "fp_common.i"

%{
#include "rayleigh_bodhaine.h"
#include "temperature.h"
%}

%base_import(rayleigh_imp_base)
%import "pressure.i"
%import "altitude.i"
%import "constant.i"

%fp_shared_ptr(FullPhysics::RayleighBodhaine);

namespace FullPhysics {

class RayleighBodhaine: public RayleighImpBase {
public:
    RayleighBodhaine(const boost::shared_ptr<Pressure>& Pres,
                     const std::vector<boost::shared_ptr<Altitude> >& Alt,
                     const boost::shared_ptr<Constant>& C);

    virtual DoubleWithUnit cross_section(const DoubleWithUnit& W) const;

    virtual boost::shared_ptr<Rayleigh> clone() const;

    virtual void print(std::ostream& Os) const;

private:
};
}
