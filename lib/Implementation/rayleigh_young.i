%include "fp_common.i"

%{
#include "rayleigh_young.h"
#include "temperature.h"
%}

%base_import(rayleigh_imp_base)

%fp_shared_ptr(FullPhysics::RayleighYoung);

namespace FullPhysics {

class RayleighYoung: public RayleighImpBase {
public:
    RayleighYoung(const boost::shared_ptr<Pressure>& Pres,
                  const std::vector<boost::shared_ptr<Altitude> >& Alt,
                  const boost::shared_ptr<Constant>& C);

    virtual DoubleWithUnit cross_section(const DoubleWithUnit& W) const;

    virtual boost::shared_ptr<Rayleigh> clone() const;

    virtual void print(std::ostream& Os) const;

private:
};
}
