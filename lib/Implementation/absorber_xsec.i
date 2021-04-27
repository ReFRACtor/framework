%include "fp_common.i"

%{
#include "absorber_xsec.h"
#include "state_mapping_linear.h"
#include "sub_state_vector_observer.h"
%}

%base_import(absorber)

%import "xsec_table.i"

%fp_shared_ptr(FullPhysics::AbsorberXSec);

namespace FullPhysics {

class AbsorberXSec: virtual public Absorber {
 public:
    AbsorberXSec(const std::vector<boost::shared_ptr<AbsorberVmr> > Vmr,
                 const boost::shared_ptr<Pressure>& Press,
                 const boost::shared_ptr<Temperature>& Temp,
                 const std::vector<boost::shared_ptr<Altitude> >& Alt,
                 const std::vector<boost::shared_ptr<XSecTable> >& XSec_tables);

    virtual ~AbsorberXSec() = default;

    virtual int number_species() const;
    virtual int number_spectrometer() const;
    virtual int number_layer() const;
    virtual std::string gas_name(int Species_index) const;

    %python_attribute(total_air_number_density_level, ArrayAdWithUnit<double, 1>)
    %python_attribute(gas_number_density_level, ArrayAdWithUnit<double, 2>)

    virtual ArrayAdWithUnit<double, 1> total_air_number_density_layer(int spec_index) const;
    virtual ArrayAdWithUnit<double, 2> gas_number_density_layer(int spec_index) const;

    virtual ArrayAd<double, 2> optical_depth_each_layer(double wn, int spec_index) const;

    virtual boost::shared_ptr<Absorber> clone() const;
    virtual boost::shared_ptr<AbsorberVmr> absorber_vmr(const std::string& Gas_name) const;
    virtual boost::shared_ptr<XSecTable> xsec_table(const std::string& Gas_name) const;

    const Pressure& pressure() const;

    %pickle_serialization();
};
}
