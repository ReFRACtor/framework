%include "fp_common.i"

%{
#include "oss_forward_model.h"
%}
%base_import(forward_model)
%base_import(named_spectrum)
%import "absorber_vmr.i"
%import "temperature.i"
%import "surface_temperature.i"
%import "ground_piecewise.i"
%import "spectral_domain.i"
%import "spectrum.i"
%import "oss_interface.i"
%import "oss_retrieval_flags.i"
%import "observer.i"


%fp_shared_ptr(FullPhysics::OssForwardModel);
%fp_shared_ptr(FullPhysics::OssFixedInputs);
%fp_shared_ptr(FullPhysics::OssMasters);
%fp_shared_ptr(FullPhysics::OssModifiedOutputs);

namespace FullPhysics {

class OssForwardModel: public ForwardModel, public Observable<boost::shared_ptr<NamedSpectrum> >  {
public:
    OssForwardModel(std::vector<boost::shared_ptr<AbsorberVmr>>& Vmr,
            const boost::shared_ptr<Pressure>& Pressure_,
            const boost::shared_ptr<Temperature>& Temperature_,
            const boost::shared_ptr<SurfaceTemperature>& Skin_temperature,
            const boost::shared_ptr<GroundPiecewise>& Ground_,
            DoubleWithUnit Obs_zen_ang, DoubleWithUnit Sol_zen_ang,
            DoubleWithUnit Lat, DoubleWithUnit Surf_alt, bool Lambertian,
            const std::string& Sel_file, const std::string& Od_file, const std::string& Sol_file,
            const std::string& Fix_file, const std::string& Ch_sel_file,
            std::vector<boost::shared_ptr<SpectralDomain>> channel_domains =
                    std::vector<boost::shared_ptr<SpectralDomain>>(),
            int Max_chans = 20000);
    virtual ~OssForwardModel();
    virtual void setup_grid();
    virtual int num_channels() const;
    virtual SpectralDomain spectral_domain(int Spec_index) const;
    virtual SpectralDomain::TypePreference spectral_domain_type_preference() const;
    virtual Spectrum radiance(int channel_index, bool skip_jacobian = false) const;
    virtual void setup_retrieval(const boost::shared_ptr<OssRetrievalFlags>& Retrieval_flags);
    virtual void print(std::ostream& Os) const;
    virtual void add_observer(Observer<boost::shared_ptr<NamedSpectrum> > & Obs);
    virtual void remove_observer(Observer<boost::shared_ptr<NamedSpectrum> >& Obs);
    void notify_spectrum_update(const Spectrum& updated_spec, const std::string& spec_name, int channel_index) const;    
    boost::shared_ptr<OssFixedInputs> fixed_inputs;
    boost::shared_ptr<OssMasters> oss_master;
    mutable boost::shared_ptr<OssModifiedOutputs> cached_outputs;
    mutable boost::shared_ptr<OssRetrievalFlags> retrieval_flags;
};
}
    
