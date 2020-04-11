#ifndef OSS_INTERFACE_H
#define OSS_INTERFACE_H

#include <string>
#include <vector>
#include <blitz/array.h>

#include "generic_object.h"
#include "array_with_unit.h"
#include "float_with_unit.h"

extern "C" {
void cppinitwrapper(int &, int &, char *, int &, int &, char *, const char *,
        int &, const char *, int &, const char *, int &, const char *, int &,
        const char *, int &, int &, int &, float &, int &, float *,
        const int &);

void cppfwdwrapper(int &, int &, float *, float *, float &, float *, int &,
        float *, float *, float &, float &, int &, float *, float *, float *,
        float &, float &, float &, float &, int &, int &, int &, float *,
        float *, float *, float *, float *, float *, float *, float *);
}

namespace FullPhysics {

/****************************************************************//**
 This class helps manage OSS 1d C strings created from components in
 a vector<std::string>
 *******************************************************************/
class OssString: public virtual GenericObject {
public:
	OssString() { } /* Allow uninitialized object creation */
	OssString(const std::vector<std::string>& Components);
	int substrlen; ///< Length of individual strings (with padding) within larger concat'd string
	int num_substr; ///< Number of individual strings within larger concat'd string
	std::string oss_str; ///< OSS specific concat'd string with padding to max substrs fixed length
	std::vector<std::string> components; ///< Backing container of OSS strings

private:
    int max_substrlen(std::vector<std::string> names) const;
    std::string str_vec_to_oss_str(std::vector<std::string> names) const;
};

/****************************************************************//**
 This class stores the dynamic outputs set by OSS each FM call
 *******************************************************************/
class OssModifiedOutputs: public virtual GenericObject {
public:
    OssModifiedOutputs() { } /* Allow creation of uninitialized objects */
    OssModifiedOutputs(blitz::Array<float, 1>& Y, blitz::Array<float, 1>& Xk_temp, blitz::Array<float, 1>& Xk_tskin,
            blitz::Array<float, 1>& Xk_out_gas, blitz::Array<float, 1>& Xk_em, blitz::Array<float, 1>& Xk_rf,
            blitz::Array<float, 1>& Xk_cldln_pres, blitz::Array<float, 1>& Xk_cldln_ext) :
            y(Y, Unit("W m^{-2} sr^{-1} cm^{-1}")), xk_temp(Xk_temp, Unit("W m^{-2} sr^{-1} cm^{-1} K^{-1}")),
			xk_tskin(Xk_tskin, Unit("W m^{-2} sr^{-1} cm^{-1} K^{-1}")),
			xk_out_gas(Xk_out_gas, Unit("W m^{-2} sr^{-1} cm^{-1}")),
			xk_em(Xk_em, Unit("W / m^2 / sr / cm")), xk_rf(Xk_rf,Unit("W m^{-2} sr^{-1} cm^{-1}")),
			xk_cldln_pres(Xk_cldln_pres, Unit("W m^{-2} sr^{-1} cm^{-1}")),
			xk_cldln_ext(Xk_cldln_ext, Unit("W m^{-2} sr^{-1} cm^{-1}")){
    }

    OssModifiedOutputs(ArrayWithUnit<float, 1>& Y, ArrayWithUnit<float, 1>& Xk_temp, ArrayWithUnit<float, 1>& Xk_tskin,
            ArrayWithUnit<float, 1>& Xk_out_gas, ArrayWithUnit<float, 1>& Xk_em, ArrayWithUnit<float, 1>& Xk_rf,
            ArrayWithUnit<float, 1>& Xk_cldln_pres, ArrayWithUnit<float, 1>& Xk_cldln_ext) :
            y(Y), xk_temp(Xk_temp), xk_tskin(Xk_tskin), xk_out_gas(Xk_out_gas), xk_em(Xk_em),
            xk_rf(Xk_rf), xk_cldln_pres(Xk_cldln_pres), xk_cldln_ext(Xk_cldln_ext){
    }

    ArrayWithUnit<float, 1> y; ///< radiance (W m−2 str−1 cm−1)
    ArrayWithUnit<float, 1> xk_temp; ///< temperature profile Jacobians (W m−2 str−1 cm−1 K−1)
    ArrayWithUnit<float, 1> xk_tskin; ///< skin temperature Jacobians (W m−2 str−1 cm−1 K−1)
    ArrayWithUnit<float, 1> xk_out_gas; ///< gas log concentration Jacobians (W m−2 str−1 cm−1)
    ArrayWithUnit<float, 1> xk_em; ///< emissivity Jacobians (W m−2 str−1 cm−1)
    ArrayWithUnit<float, 1> xk_rf; ///< reflectivity Jacobians (W m−2 str−1 cm−1)
    ArrayWithUnit<float, 1> xk_cldln_pres; ///< cloud center log pressure Jacobians  (W m−2 str−1 cm−1)
    ArrayWithUnit<float, 1> xk_cldln_ext; ///< cloud peak log extinction Jacobians  (W m−2 str−1 cm−1)
};


/****************************************************************//**
 This class stores the fixed outputs set by OSS FM on init
 *******************************************************************/
class OssFixedOutputs: public virtual GenericObject {
public:
    OssFixedOutputs() { } /* Allow uninitialized object creation */
    OssFixedOutputs(int Num_chan, blitz::Array<float, 1>& Center_wavenumber) :
            num_chan(Num_chan), center_wavenumber(Center_wavenumber, Unit("Wavenumbers")) {
    }

    OssFixedOutputs(int Num_chan, ArrayWithUnit<float, 1>& Center_wavenumber) :
            num_chan(Num_chan), center_wavenumber(Center_wavenumber) {
    }

    int num_chan; ///< Number of channels available in OSS RTM
    ArrayWithUnit<float, 1> center_wavenumber; //< Center wavenumbers of channels
};

/****************************************************************//**
 This class stores the fixed inputs used by the OSS FM (set on init)
 *******************************************************************/
class OssFixedInputs: public virtual GenericObject {
public:
    OssFixedInputs() { } /* Allow creation of uninitialized objects */
    OssFixedInputs(std::vector<std::string>& Gas_names,
            std::vector<std::string>& Gas_jacobian_names, std::string& Sel_file,
            std::string& Od_file, std::string& Sol_file, std::string& Fix_file,
            std::string& Ch_sel_file, int Num_vert_lev, int Num_surf_points,
            float Min_extinct_cld, int Max_chans = 20000) :
            gas_names(Gas_names), gas_jacobian_names(Gas_jacobian_names),
			sel_file(Sel_file), od_file(Od_file), sol_file(Sol_file), fix_file(Fix_file),
			ch_sel_file(Ch_sel_file),num_vert_lev(Num_vert_lev), num_surf_points(Num_surf_points),
			min_extinct_cld(Min_extinct_cld, Unit("km^-1")), max_chans(Max_chans), oss_gas_names(Gas_names),
			oss_gas_jacobian_names(Gas_jacobian_names) {
    }

    int max_chans; ///< Maximum number of channels
    std::vector<std::string> gas_names; ///< Molecular gas names
    OssString oss_gas_names; ///< OSS 1d str representation of gas names
    std::vector<std::string> gas_jacobian_names; ///< Molecular gas names for Jacobians
    OssString oss_gas_jacobian_names; ///< OSS 1d str representation of gas names for Jacobians
    std::string sel_file; ///< File name of OSS nodes and weights
    std::string od_file; ///< File name of optical property lookup table
    std::string sol_file; ///< File name of solar radiance
    std::string fix_file; ///< File name of default profiles for variable gases
    std::string ch_sel_file; ///< File name of list(s) of channel subsets (or "NULL")
    int num_vert_lev; ///< Number of vertical levels of state vector
    int num_surf_points; ///< Number of surface grid points
    FloatWithUnit min_extinct_cld; ///< Threshold of extinction for including cloud in RT (km−1)
};

/****************************************************************//**
 This class stores the variable inputs used by the OSS FM.
 Unlike OssFixedInputs, these inputs may vary by forward model call.
 *******************************************************************/
class OssModifiedInputs: public virtual GenericObject {
public:
    OssModifiedInputs(blitz::Array<float, 1>& Pressure,
            blitz::Array<float, 1>& Temp, float Skin_temp,
            blitz::Array<float, 2>& Vmr_gas, blitz::Array<float, 1>& Emis,
            blitz::Array<float, 1>& Refl, float Scale_cld, float Pressure_cld,
            blitz::Array<float, 1>& Ext_cld, blitz::Array<float, 1>& Surf_grid,
            blitz::Array<float, 1>& Cld_grid, float Obs_zen_ang, float Sol_zen_ang,
            float Lat, float Surf_alt, bool Lambertian) :
            pressure(Pressure, units::mbar), temp(Temp, units::K), skin_temp(Skin_temp, units::K),
			vmr_gas(Vmr_gas, units::dimensionless), emis(Emis, units::dimensionless),
			refl(Refl, units::dimensionless), scale_cld(Scale_cld, units::dimensionless),
			pressure_cld(Pressure_cld, units::mbar), ext_cld(Ext_cld, Unit("km^-1")),
			surf_grid(Surf_grid, units::inv_cm), cld_grid(Cld_grid, units::inv_cm),
			obs_zen_ang(Obs_zen_ang, units::deg), sol_zen_ang(Sol_zen_ang, units::deg),
			lat(Lat, units::deg), surf_alt(Surf_alt, units::m), lambertian(Lambertian) {
    }

    OssModifiedInputs(ArrayWithUnit<float, 1>& Pressure,
            ArrayWithUnit<float, 1>& Temp, FloatWithUnit Skin_temp,
            ArrayWithUnit<float, 2>& Vmr_gas, ArrayWithUnit<float, 1>& Emis,
            ArrayWithUnit<float, 1>& Refl, FloatWithUnit Scale_cld, FloatWithUnit Pressure_cld,
            ArrayWithUnit<float, 1>& Ext_cld, ArrayWithUnit<float, 1>& Surf_grid,
            ArrayWithUnit<float, 1>& Cld_grid, FloatWithUnit Obs_zen_ang, FloatWithUnit Sol_zen_ang,
			FloatWithUnit Lat, FloatWithUnit Surf_alt, bool Lambertian) :
            pressure(Pressure), temp(Temp), skin_temp(Skin_temp), vmr_gas(Vmr_gas),
			emis(Emis), refl(Refl), scale_cld(Scale_cld), pressure_cld(Pressure_cld),
			ext_cld(Ext_cld), surf_grid(Surf_grid), cld_grid(Cld_grid), obs_zen_ang(Obs_zen_ang),
			sol_zen_ang(Sol_zen_ang), lat(Lat), surf_alt(Surf_alt), lambertian(Lambertian) {
    }

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

/****************************************************************//**
 This class wraps the OSS FM exposed by the C wrapper
 *******************************************************************/
class OssMasters: public virtual GenericObject {
public:
    OssMasters() { } /* Allow creation of uninitialized objects */
    OssMasters(OssFixedInputs& Fixed_inputs) : fixed_inputs(Fixed_inputs) {}

    void init();
    OssModifiedOutputs run_fwd_model(OssModifiedInputs& Modified_inputs);

    OssFixedInputs fixed_inputs;
    OssFixedOutputs fixed_outputs;
};
}
#endif
