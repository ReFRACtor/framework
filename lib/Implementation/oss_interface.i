%include "fp_common.i"

%{
#include "oss_interface.h"
#include "spectral_domain.h"
%}
%import "spectral_domain.i"
%import "array_with_unit.i"
%import "float_with_unit.i"
%import "string_vector_to_char.i"

%fp_shared_ptr(FullPhysics::OssModifiedOutputs);
%fp_shared_ptr(FullPhysics::OssFixedOutputs);
%fp_shared_ptr(FullPhysics::OssFixedInputs);
%fp_shared_ptr(FullPhysics::OssModifiedInputs);
%fp_shared_ptr(FullPhysics::OssMasters);

namespace FullPhysics {

class OssModifiedOutputs: public virtual GenericObject {
public:
    OssModifiedOutputs(blitz::Array<float, 1>& Y, blitz::Array<float, 2>& Xk_temp, blitz::Array<float, 1>& Xk_tskin,
            blitz::Array<float, 3>& Xk_out_gas, blitz::Array<float, 2>& Xk_em, blitz::Array<float, 2>& Xk_rf,
            blitz::Array<float, 1>& Xk_cldln_pres, blitz::Array<float, 2>& Xk_cldln_ext);

    OssModifiedOutputs(ArrayWithUnit<float, 1>& Y, ArrayWithUnit<float, 2>& Xk_temp, ArrayWithUnit<float, 1>& Xk_tskin,
            ArrayWithUnit<float, 3>& Xk_out_gas, ArrayWithUnit<float, 2>& Xk_em, ArrayWithUnit<float, 2>& Xk_rf,
            ArrayWithUnit<float, 1>& Xk_cldln_pres, ArrayWithUnit<float, 2>& Xk_cldln_ext) :
            y(Y), xk_temp(Xk_temp), xk_tskin(Xk_tskin), xk_out_gas(Xk_out_gas), xk_em(Xk_em),
            xk_rf(Xk_rf), xk_cldln_pres(Xk_cldln_pres), xk_cldln_ext(Xk_cldln_ext);

    ArrayWithUnit<float, 1> y; ///< radiance (W m−2 str−1 cm−1)
    ArrayWithUnit<float, 2> xk_temp; ///< temperature profile Jacobians (W m−2 str−1 cm−1 K−1)
    ArrayWithUnit<float, 1> xk_tskin; ///< skin temperature Jacobians (W m−2 str−1 cm−1 K−1)
    ArrayWithUnit<float, 3> xk_out_gas; ///< gas log concentration Jacobians (W m−2 str−1 cm−1)
    ArrayWithUnit<float, 2> xk_em; ///< emissivity Jacobians (W m−2 str−1 cm−1)
    ArrayWithUnit<float, 2> xk_rf; ///< reflectivity Jacobians (W m−2 str−1 cm−1)
    ArrayWithUnit<float, 1> xk_cldln_pres; ///< cloud center log pressure Jacobians  (W m−2 str−1 cm−1)
    ArrayWithUnit<float, 2> xk_cldln_ext; ///< cloud peak log extinction Jacobians  (W m−2 str−1 cm−1)
};

class OssFixedOutputs: public virtual GenericObject {
public:
    OssFixedOutputs(int Num_chan, blitz::Array<float, 1>& Full_spectral_point);
    OssFixedOutputs(int Num_chan, ArrayWithUnit<float, 1>& Full_spectral_point);
    int num_chan; ///< Number of channels available in OSS RTM
    ArrayWithUnit<float, 1> full_spectral_point; //< Center spectral point of channels (cm−1)
    std::vector<ArrayWithUnit<float, 1>> sensor_spectral_point; //< Per sensor spectral points  (cm−1)
};

class OssFixedInputs: public virtual GenericObject {
public:
    OssFixedInputs(std::vector<std::string>& Gas_names,
            std::vector<std::string>& Gas_jacobian_names, std::string& Sel_file,
            std::string& Od_file, std::string& Sol_file, std::string& Fix_file,
            std::string& Ch_sel_file, int Num_vert_lev, int Num_surf_points,
            float Min_extinct_cld,
            std::vector<boost::shared_ptr<SpectralDomain>> Channel_domains =
                    std::vector<boost::shared_ptr<SpectralDomain>>(),
            int Max_chans = 20000);

    std::vector<std::string> gas_names; ///< Molecular gas names
    boost::shared_ptr<StringVectorToChar> oss_gas_names; ///< OSS 1d str representation of gas names
    std::vector<std::string> gas_jacobian_names; ///< Molecular gas names for Jacobians
    boost::shared_ptr<StringVectorToChar> oss_gas_jacobian_names; ///< OSS 1d str representation of gas names for Jacobians
    std::string sel_file; ///< File name of OSS nodes and weights
    std::string od_file; ///< File name of optical property lookup table
    std::string sol_file; ///< File name of solar radiance
    std::string fix_file; ///< File name of default profiles for variable gases
    std::string ch_sel_file; ///< File name of list(s) of channel subsets (or "NULL")
    int num_vert_lev; ///< Number of vertical levels of state vector
    int num_surf_points; ///< Number of surface grid points
    FloatWithUnit min_extinct_cld; ///< Threshold of extinction for including cloud in RT (km−1)
	std::vector<boost::shared_ptr<SpectralDomain>> channel_domains;
    int max_chans; ///< Maximum number of channels
};

class OssModifiedInputs: public virtual GenericObject {
public:
    OssModifiedInputs(blitz::Array<float, 1>& Pressure,
            blitz::Array<float, 1>& Temp, float Skin_temp,
            blitz::Array<float, 2>& Vmr_gas, blitz::Array<float, 1>& Emis,
            blitz::Array<float, 1>& Refl, float Scale_cld, float Pressure_cld,
            blitz::Array<float, 1>& Ext_cld, blitz::Array<float, 1>& Surf_grid,
            blitz::Array<float, 1>& Cld_grid, float Obs_zen_ang, float Sol_zen_ang,
            float Lat, float Surf_alt, bool Lambertian);

    OssModifiedInputs(ArrayWithUnit<float, 1>& Pressure,
            ArrayWithUnit<float, 1>& Temp, FloatWithUnit Skin_temp,
            ArrayWithUnit<float, 2>& Vmr_gas, ArrayWithUnit<float, 1>& Emis,
            ArrayWithUnit<float, 1>& Refl, FloatWithUnit Scale_cld, FloatWithUnit Pressure_cld,
            ArrayWithUnit<float, 1>& Ext_cld, ArrayWithUnit<float, 1>& Surf_grid,
            ArrayWithUnit<float, 1>& Cld_grid, FloatWithUnit Obs_zen_ang, FloatWithUnit Sol_zen_ang,
			FloatWithUnit Lat, FloatWithUnit Surf_alt, bool Lambertian);

    ArrayWithUnit<float, 1> pressure; //< Atmospheric pressure profile (ordered from high to low) (mb)
    ArrayWithUnit<float, 1> temp; //< Temperature profile (K)
    FloatWithUnit skin_temp; //< Skin temperature (K)
    ArrayWithUnit<float, 2> vmr_gas; //< Volume mixing ratio of gases (unitless)
    ArrayWithUnit<float, 1> emis; //< Surface emissivity (unitless)
    ArrayWithUnit<float, 1> refl; //< Surface reflectivity (unitless)
    FloatWithUnit scale_cld; //< Scale of cloud log thickness
    FloatWithUnit pressure_cld; //< Pressure at cloud center (peak of extinction profile) (mb)
    ArrayWithUnit<float, 1> ext_cld; //< Cloud peak extinction (km−1)
    ArrayWithUnit<float, 1> surf_grid; //< Vector of surface grid point wavenumbers (cm−1)
    ArrayWithUnit<float, 1> cld_grid; //< Cloud grid point wavenumbers (cm−1)
    FloatWithUnit obs_zen_ang; //< Observation zenith angle (degrees)
    FloatWithUnit sol_zen_ang; //< Solar zenith angle (degrees)
    FloatWithUnit lat; //< latitude (degrees)
    FloatWithUnit surf_alt; //< surface altitude (m)
    bool lambertian; //< flag for Lambertian (true) or (default,false) specular surface*/
};

class OssMasters: public virtual GenericObject {
public:
    OssMasters(boost::shared_ptr<OssFixedInputs> Fixed_inputs);
    ~OssMasters();
    void init(double spec_domain_tolerance=0.001);
    boost::shared_ptr<OssModifiedOutputs> run_fwd_model(int Channel_index,
        boost::shared_ptr<OssModifiedInputs> Modified_inputs);

    boost::shared_ptr<OssFixedInputs> fixed_inputs;
    boost::shared_ptr<OssFixedOutputs> fixed_outputs;
};
}

%template(vector_spectral_domain) std::vector<boost::shared_ptr<FullPhysics::SpectralDomain> >;
