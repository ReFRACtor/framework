#ifndef OSS_FORWARD_MODEL_H
#define OSS_FORWARD_MODEL_H

#include "forward_model.h"
#include "absorber_vmr.h"
#include "temperature.h"
#include "surface_temperature.h"
#include "ground_piecewise.h"
#include "oss_interface.h"
#include "oss_retrieval_flags.h"

#include <string>

namespace FullPhysics {
/****************************************************************//**
  This a forward model class that wraps the AER OSS Forward Model
*******************************************************************/

class OssForwardModel : public ForwardModel {
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
    virtual ~OssForwardModel() {}
    virtual void setup_grid();

    virtual int num_channels() const { return 1; }

    virtual SpectralDomain spectral_domain(int Spec_index) const {
        return SpectralDomain(center_spectral_point);
    }
    virtual SpectralDomain::TypePreference spectral_domain_type_preference() const {
        return SpectralDomain::PREFER_WAVENUMBER;
    }
    virtual Spectrum radiance(int channel_index, bool skip_jacobian = false) const;
    virtual void setup_retrieval(OssRetrievalFlags Retrieval_flags) const;
    virtual void print(std::ostream& Os) const { Os << "OssForwardModel"; }

    boost::shared_ptr<OssFixedInputs> fixed_inputs;
    boost::shared_ptr<OssMasters> oss_master;
    mutable boost::shared_ptr<OssModifiedOutputs> cached_outputs;
    mutable boost::shared_ptr<OssRetrievalFlags> retrieval_flags;
private:
    std::vector<boost::shared_ptr<AbsorberVmr>> vmr;
    boost::shared_ptr<Pressure> pressure;
    boost::shared_ptr<Temperature> temperature;
    boost::shared_ptr<SurfaceTemperature> skin_temperature;
    boost::shared_ptr<GroundPiecewise> ground;
    DoubleWithUnit obs_zen_ang;
    DoubleWithUnit sol_zen_ang;
    DoubleWithUnit lat;
    DoubleWithUnit surf_alt;
    bool lambertian;

    std::string sel_file;
    int sel_file_sz;
    std::string od_file;
    int od_file_sz;
    std::string sol_file;
    int sol_file_sz;
    std::string fix_file;
    int fix_file_sz;
    std::string ch_sel_file;
    int ch_sel_file_sz;
    int max_chans;

    std::vector<boost::shared_ptr<SpectralDomain>> channel_domains;
    ArrayWithUnit<double, 1> center_spectral_point;
};
}
#endif
