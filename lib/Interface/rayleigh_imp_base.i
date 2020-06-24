%include "fp_common.i"

%{
#include "rayleigh_imp_base.h"
#include "temperature.h"
%}

%base_import(pressure)
%base_import(altitude)
%base_import(rayleigh)

%import "constant.i"
%import "double_with_unit.i"

%fp_shared_ptr(FullPhysics::RayleighImpBase);
%fp_shared_ptr(FullPhysics::Observable<FullPhysics::RayleighImpBase>);
%fp_shared_ptr(FullPhysics::Observer<FullPhysics::RayleighImpBase>);

namespace FullPhysics {

class RayleighImpBase: public Observer<Pressure>, public Observer<Altitude>, public Rayleigh {
public:
    RayleighImpBase(const boost::shared_ptr<Pressure>& Pres, 
                    const std::vector<boost::shared_ptr<Altitude> >& Alt,
                    const Constant& C);
    virtual void notify_update(const Pressure& P);
    virtual void notify_update(const Altitude& A);
    virtual ArrayAd<double, 1> optical_depth_each_layer(double wn, int spec_index) const;
    virtual DoubleWithUnit cross_section(const DoubleWithUnit& W) const = 0;
    virtual void print(std::ostream& Os) const;
  %pickle_serialization();

};

%template(ObservableRayleighImpBase) FullPhysics::Observable<RayleighImpBase>;
%template(ObserverRayleighImpBase) FullPhysics::Observer<RayleighImpBase>;

}
