%include "fp_common.i"

%{
#include "oss_forward_model.h"
%}

%fp_shared_ptr(FullPhysics::OssForwardModel);

namespace FullPhysics {

class OssForwardModel: public virtual GenericObject {
public:
    OssForwardModel(std::vector<boost::shared_ptr<AbsorberVmr>>& Vmr,
            const boost::shared_ptr<Pressure>& Pressure_,
            const boost::shared_ptr<Temperature>& Temperature_,
            const boost::shared_ptr<SurfaceTemperature>& Skin_temperature,
            const boost::shared_ptr<GroundPiecewise>& Ground_,
            DoubleWithUnit Obs_zen_ang, DoubleWithUnit Sol_zen_ang,
            DoubleWithUnit Lat, DoubleWithUnit Surf_alt, bool Lambertian,
            const std::string& Sel_file, const std::string& Od_file, const std::string& Sol_file,
            const std::string& Fix_file, const std::string& Ch_sel_file, int Max_chans = 20000);
    virtual ~OssForwardModel();
    virtual void setup_grid();
    virtual int num_channels() const;
    virtual SpectralDomain spectral_domain(int Spec_index) const;
    virtual SpectralDomain::TypePreference spectral_domain_type_preference() const;
    virtual boost::shared_ptr<StateVector> state_vector() const;
    virtual Spectrum radiance(int channel_index, bool skip_jacobian = false) const;
    virtual void print(std::ostream& Os) const;
};
}
    