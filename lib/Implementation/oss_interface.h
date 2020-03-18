#ifndef OSS_INTERFACE_H
#define OSS_INTERFACE_H

#include <string>
#include <vector>
#include <blitz/array.h>

#include "generic_object.h"
#include "array_with_unit.h"

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
 This class stores the dynamic outputs set by OSS each FM call
 *******************************************************************/
class OssModifiedOutputs: public virtual GenericObject {
public:
    /* TODO: Remove default constructor */
    OssModifiedOutputs() {}
    OssModifiedOutputs(blitz::Array<float, 1>& Y, blitz::Array<float, 1>& XkTemp, blitz::Array<float, 1>& XkTskin,
            blitz::Array<float, 1>& XkOutGas, blitz::Array<float, 1>& XkEm, blitz::Array<float, 1>& XkRf,
            blitz::Array<float, 1>& XkCldlnPres, blitz::Array<float, 1>& XkCldlnExt) :
            y(Y), xk_temp(XkTemp), xk_tskin(XkTskin), xk_out_gas(XkOutGas), xk_em(XkEm),
            xk_rf(XkRf), xk_cldln_pres(XkCldlnPres), xk_cldln_ext(XkCldlnExt){
    }

    const blitz::Array<float, 1>& XkCldlnExt() const {
        return xk_cldln_ext;
    }

    const blitz::Array<float, 1>& XkCldlnPres() const {
        return xk_cldln_pres;
    }

    const blitz::Array<float, 1>& XkEm() const {
        return xk_em;
    }

    const blitz::Array<float, 1>& XkOutGas() const {
        return xk_out_gas;
    }

    const blitz::Array<float, 1>& XkRf() const {
        return xk_rf;
    }

    const blitz::Array<float, 1>& XkTemp() const {
        return xk_temp;
    }

    const blitz::Array<float, 1>& XkTskin() const {
        return xk_tskin;
    }

    const blitz::Array<float, 1>& Y() const {
        return y;
    }

private:
    /* TODO: ArrayWithUnit */
    blitz::Array<float, 1> y; ///< radiance (W m−2 str−1 cm−1)
    blitz::Array<float, 1> xk_temp; ///< temperature profile Jacobians (W m−2 str−1 cm−1 K−1)
    blitz::Array<float, 1> xk_tskin; ///< skin temperature Jacobians (W m−2 str−1 cm−1 K−1)
    blitz::Array<float, 1> xk_out_gas; ///< gas log concentration Jacobians (W m−2 str−1 cm−1)
    blitz::Array<float, 1> xk_em; ///< emissivity Jacobians (W m−2 str−1 cm−1)
    blitz::Array<float, 1> xk_rf; ///< reflectivity Jacobians (W m−2 str−1 cm−1)
    blitz::Array<float, 1> xk_cldln_pres; ///< cloud center log pressure Jacobians  (W m−2 str−1 cm−1)
    blitz::Array<float, 1> xk_cldln_ext; ///< cloud peak log extinction Jacobians  (W m−2 str−1 cm−1)
};


/****************************************************************//**
 This class stores the fixed outputs set by OSS FM on init
 *******************************************************************/
class OssFixedOutputs: public virtual GenericObject {
public:
	/* Constructor to create a empty object */
	OssFixedOutputs() {}
    OssFixedOutputs(int NumChan, ArrayWithUnit<float, 1>& CenterWavenumber) :
            num_chan(NumChan), center_wavenumber(CenterWavenumber) {
    }

    const ArrayWithUnit<float, 1>& CenterWavenumber() const {
        return center_wavenumber;
    }

    int NumChan() const {
        return num_chan;
    }

private:
    int num_chan; ///< Number of channels available in OSS RTM
    ArrayWithUnit<float, 1> center_wavenumber; //< Center wavenumbers of channels

};

/****************************************************************//**
 This class stores the fixed inputs used by the OSS FM (set on init)
 *******************************************************************/
class OssFixedInputs: public virtual GenericObject {
public:
    OssFixedInputs(int NumMolecules, std::vector<std::string>& GasNames,
            std::vector<std::string>& GasJacobianNames, std::string& SelFile,
            std::string& OdFile, std::string& SolFile, std::string& FixFile,
            std::string& ChSelFile, int NumVertLev, int NumSurfPoints,
            float MinExtinctCld, int MaxChans = 20000) :
            num_molecules(NumMolecules), gas_names(GasNames), gas_jacobian_names(
                    GasJacobianNames), sel_file(SelFile), od_file(OdFile), sol_file(
                    SolFile), fix_file(FixFile), ch_sel_file(ChSelFile), num_vert_lev(
                    NumVertLev), num_surf_points(NumSurfPoints), min_extinct_cld(
                    MinExtinctCld), max_chans(MaxChans) {

    }

    virtual ~OssFixedInputs() = default;

    const std::string& ChSelFile() const {
        return ch_sel_file;
    }

    const std::string& FixFile() const {
        return fix_file;
    }

    const std::vector<std::string>& GasJacobianNames() const {
        /* TODO: Translate to 1d C char array */
        return gas_jacobian_names;
    }

    const std::vector<std::string>& GasNames() const {
        /* TODO: Translate to 1d char array */
        return gas_names;
    }

    int MaxChans() const {
        return max_chans;
    }

    float MinExtinctCld() const {
        return min_extinct_cld;
    }

    int NumMolecules() const {
        return num_molecules;
    }

    int NumSurfPoints() const {
        return num_surf_points;
    }

    int NumVertLev() const {
        return num_vert_lev;
    }

    const std::string& OdFile() const {
        return od_file;
    }

    const std::string& SelFile() const {
        return sel_file;
    }

    const std::string& SolFile() const {
        return sol_file;
    }

private:
    int max_chans; ///< Maximum number of channels
    int num_molecules; ///< Number of molecules in input profiles
    std::vector<std::string> gas_names; ///< Molecular gas names
    std::vector<std::string> gas_jacobian_names; ///< Molecular gas names for Jacobians
    std::string sel_file; ///< File name of OSS nodes and weights
    std::string od_file; ///< File name of optical property lookup table
    std::string sol_file; ///< File name of solar radiance
    std::string fix_file; ///< File name of default profiles for variable gases
    std::string ch_sel_file; ///< File name of list(s) of channel subsets (or "NULL")
    int num_vert_lev; ///< Number of vertical levels of state vector
    int num_surf_points; ///< Number of surface grid points
    float min_extinct_cld; ///< Threshold of extinction for including cloud in RT (km−1)
};

/****************************************************************//**
 This class stores the variable inputs used by the OSS FM.
 Unlike OssFixedInputs, these inputs may vary by forward model call.
 *******************************************************************/
class OssModifiedInputs: public virtual GenericObject {
public:
    OssModifiedInputs(blitz::Array<float, 1>& Pressure,
            blitz::Array<float, 1>& Temp, float SkinTemp,
            blitz::Array<float, 2>& VmrGas, blitz::Array<float, 1>& Emis,
            blitz::Array<float, 1>& Refl, float ScaleCld, float PressureCld,
            blitz::Array<float, 1>& ExtCld, blitz::Array<float, 1>& SurfGrid,
            blitz::Array<float, 1>& CldGrid, float ObsZenAng, float SolZenAng,
            float Lat, float SurfAlt, bool Lambertian) :
            pressure(Pressure), temp(Temp), skin_temp(SkinTemp), vmr_gas(
                    VmrGas), emis(Emis), refl(Refl), scale_cld(ScaleCld), pressure_cld(
                    PressureCld), ext_cld(ExtCld), surf_grid(SurfGrid), cld_grid(
                    CldGrid), obs_zen_ang(ObsZenAng), sol_zen_ang(SolZenAng), lat(
                    Lat), surf_alt(SurfAlt), lambertian(Lambertian) {
    }

    virtual ~OssModifiedInputs() = default;

    blitz::Array<float, 1>& CldGrid() {
        return cld_grid;
    }

    void setCldGrid(const blitz::Array<float, 1>& cldGrid) {
        cld_grid = cldGrid;
    }

    blitz::Array<float, 1>& Emis() {
        return emis;
    }

    void setEmis(const blitz::Array<float, 1>& emis) {
        this->emis = emis;
    }

    blitz::Array<float, 1>& ExtCld() {
        return ext_cld;
    }

    void setExtCld(const blitz::Array<float, 1>& extCld) {
        ext_cld = extCld;
    }

    bool isLambertian() const {
        return lambertian;
    }

    void setLambertian(bool lambertian) {
        this->lambertian = lambertian;
    }

    float Lat() const {
        return lat;
    }

    void setLat(float lat) {
        this->lat = lat;
    }

    float ObsZenAng() const {
        return obs_zen_ang;
    }

    void setObsZenAng(float obsZenAng) {
        obs_zen_ang = obsZenAng;
    }

    blitz::Array<float, 1>& Pressure() {
        return pressure;
    }

    void setPressure(const blitz::Array<float, 1>& pressure) {
        this->pressure = pressure;
    }

    float PressureCld() const {
        return pressure_cld;
    }

    void setPressureCld(float pressureCld) {
        pressure_cld = pressureCld;
    }

    blitz::Array<float, 1>& Refl() {
        return refl;
    }

    void setRefl(const blitz::Array<float, 1>& refl) {
        this->refl = refl;
    }

    float ScaleCld() const {
        return scale_cld;
    }

    void setScaleCld(float scaleCld) {
        scale_cld = scaleCld;
    }

    float SkinTemp() const {
        return skin_temp;
    }

    void setSkinTemp(float skinTemp) {
        skin_temp = skinTemp;
    }

    float SolZenAng() const {
        return sol_zen_ang;
    }

    void setSolZenAng(float solZenAng) {
        sol_zen_ang = solZenAng;
    }

    float SurfAlt() const {
        return surf_alt;
    }

    void setSurfAlt(float surfAlt) {
        surf_alt = surfAlt;
    }

    blitz::Array<float, 1>& SurfGrid() {
        return surf_grid;
    }

    void setSurfGrid(const blitz::Array<float, 1>& surfGrid) {
        surf_grid = surfGrid;
    }

    blitz::Array<float, 1>& Temp() {
        return temp;
    }

    void setTemp(const blitz::Array<float, 1>& temp) {
        this->temp = temp;
    }

    blitz::Array<float, 2>& VmrGas() {
        return vmr_gas;
    }

    void setVmrGas(const blitz::Array<float, 2>& vmrGas) {
        vmr_gas = vmrGas;
    }

private:
    /* TODO: ArrayWithUnit and convert on use instead? */
    blitz::Array<float, 1> pressure; //< Atmospheric pressure profile (ordered from high to low) (mb)
    blitz::Array<float, 1> temp; //< Temperature profile (K)
    float skin_temp; //< Skin temperature (K)
    blitz::Array<float, 2> vmr_gas; //< Volume mixing ratio of gases (unitless)
    blitz::Array<float, 1> emis; //< Surface emissivity (unitless)
    blitz::Array<float, 1> refl; //< Surface reflectivity (unitless)
    float scale_cld; //< Scale of cloud log thickness
    float pressure_cld; //< Pressure at cloud center (peak of extinction profile) (mb)
    blitz::Array<float, 1> ext_cld; //< Cloud peak extinction (km−1)
    blitz::Array<float, 1> surf_grid; //< Vector of surface grid point wavenumbers (cm−1)
    blitz::Array<float, 1> cld_grid; //< Cloud grid point wavenumbers (cm−1)
    float obs_zen_ang; //< Observation zenith angle (degrees)
    float sol_zen_ang; //< Solar zenith angle (degrees)
    float lat; //< latitude (degrees)
    float surf_alt; //< surface altitude (m)
    bool lambertian; //< flag for Lambertian (true) or (default,false) specular surface*/
};

/* OssMasters (boost::shared_ptr<OssFixedInputs>)
 *
 *
 * init() -> wraps cppinitwrapper (sets internal OssInitOutputs)
 *
 * run_fwd_model(boost::shared_ptr<OssModifiedInputs>) -> wraps cppfwdwrapper (returns OssFmOutputs)
 *
 * boost::shared_ptr<OssInitOutputs> getInitInputs()
 *
 *
 * */

/****************************************************************//**
 This class wraps the OSS FM exposed by the C wrapper
 *******************************************************************/
class OssMasters: public virtual GenericObject {
public:
    OssMasters(OssFixedInputs& FixedInputs) :
        fixed_inputs(FixedInputs) {}

    void init() {
        ArrayWithUnit<float, 1> center_wavenumbers = ArrayWithUnit<float, 1>(blitz::Array<float, 1>(fixed_inputs.MaxChans()),
                units::inv_cm);

        int num_chan;
        // create OSS 1d gas string from vec
        int num_gases = fixed_inputs.GasNames().size();
        int len_gas_name = max_strlen(fixed_inputs.GasNames());
        std::string oss_gas_names = str_vec_to_oss_str(fixed_inputs.GasNames());
        // create OSS 1d gas jacobian string from vec
        int num_jacob = fixed_inputs.GasJacobianNames().size();
        int len_jacob_name = max_strlen(fixed_inputs.GasJacobianNames());
        std::string oss_jacob_gas_names = str_vec_to_oss_str(fixed_inputs.GasJacobianNames());

        // better way?
        int sel_file_sz = fixed_inputs.SelFile().length();
        int od_file_sz  = fixed_inputs.OdFile().length();
        int sol_file_sz = fixed_inputs.SolFile().length();
        int fix_file_sz = fixed_inputs.FixFile().length();
        int ch_sel_file_sz = fixed_inputs.ChSelFile().length();
        int num_vert_lev = fixed_inputs.NumVertLev();
        int num_surf_points = fixed_inputs.NumSurfPoints();
        float min_extinct_cld = fixed_inputs.MinExtinctCld();

        cppinitwrapper(num_gases, len_gas_name, &oss_gas_names[0],
                num_jacob, len_jacob_name, &oss_jacob_gas_names[0],
                fixed_inputs.SelFile().c_str(), sel_file_sz,
                fixed_inputs.OdFile().c_str(), od_file_sz,
                fixed_inputs.SolFile().c_str(), sol_file_sz,
                fixed_inputs.FixFile().c_str(), fix_file_sz,
                fixed_inputs.ChSelFile().c_str(), ch_sel_file_sz,
                num_vert_lev, num_surf_points,
                min_extinct_cld, num_chan,
                center_wavenumbers.value.data(), fixed_inputs.MaxChans());

        center_wavenumbers.value.resizeAndPreserve(num_chan);

        fixed_outputs = OssFixedOutputs(num_chan, center_wavenumbers);
    }

    OssModifiedOutputs run_fwd_model(OssModifiedInputs& ModifiedInputs) {
        /* Outputs */
        blitz::Array<float, 1> y = blitz::Array<float, 1>(fixed_outputs.NumChan());
        blitz::Array<float, 1> xk_temp = blitz::Array<float, 1>(fixed_outputs.NumChan() * fixed_inputs.NumVertLev());
        blitz::Array<float, 1> xk_tskin = blitz::Array<float, 1>(fixed_outputs.NumChan());
        blitz::Array<float, 1> xk_out_gas = blitz::Array<float, 1>(fixed_outputs.NumChan() * fixed_inputs.NumVertLev() *
                fixed_inputs.GasJacobianNames().size());
        blitz::Array<float, 1> xk_em = blitz::Array<float, 1>(fixed_outputs.NumChan() * fixed_inputs.NumSurfPoints());
        blitz::Array<float, 1> xk_rf = blitz::Array<float, 1>(fixed_outputs.NumChan() * fixed_inputs.NumSurfPoints());
        blitz::Array<float, 1> xk_cldln_pres = blitz::Array<float, 1>(fixed_outputs.NumChan());
        blitz::Array<float, 1> xk_cldln_ext = blitz::Array<float, 1>(fixed_outputs.NumChan() * ModifiedInputs.CldGrid().rows());

        /* Inputs to use as l-values. Avoids: "expects an l-value for ... argument" */
        int num_vert_lev = fixed_inputs.NumVertLev();
        int num_gas = fixed_inputs.GasNames().size();
        float skin_temp = ModifiedInputs.SkinTemp();
        int num_surf_grid = ModifiedInputs.SurfGrid().rows();
        float scale_cld = ModifiedInputs.ScaleCld();
        float pressure_cld = ModifiedInputs.PressureCld();
        int num_cld = ModifiedInputs.CldGrid().rows();
        float obs_zen_ang = ModifiedInputs.ObsZenAng();
        float sol_zen_ang = ModifiedInputs.SolZenAng();
        float lat = ModifiedInputs.Lat();
        float surf_alt = ModifiedInputs.SurfAlt();
        int is_lambertian = ModifiedInputs.isLambertian();
        int num_gas_jacob = fixed_inputs.GasJacobianNames().size();
        int num_chan = fixed_outputs.NumChan();

        cppfwdwrapper(num_vert_lev, num_gas, ModifiedInputs.Pressure().data(),
        		ModifiedInputs.Temp().data(), skin_temp, ModifiedInputs.VmrGas().data(),
				num_surf_grid, ModifiedInputs.Emis().data(), ModifiedInputs.Refl().data(),
				scale_cld, pressure_cld, num_cld,
                ModifiedInputs.ExtCld().data(), ModifiedInputs.SurfGrid().data(), ModifiedInputs.CldGrid().data(),
				obs_zen_ang, sol_zen_ang, lat,
				surf_alt, is_lambertian, num_gas_jacob,
				num_chan, y.data(), xk_temp.data(),
				xk_tskin.data(), xk_out_gas.data(), xk_em.data(),
				xk_rf.data(), xk_cldln_pres.data(), xk_cldln_ext.data());

        return OssModifiedOutputs(y, xk_temp, xk_tskin, xk_out_gas, xk_em, xk_rf, xk_cldln_pres, xk_cldln_ext);
    }

    const OssFixedOutputs& FixedOutputs() const {
        return fixed_outputs;
    }


private:
    int max_strlen(std::vector<std::string> names) {
        int max_len = -1;
        for (const auto& name : names) {
            if (name.length() > max_len) {
                max_len = name.length();
            }
        }
        return max_len;
    }

    std::string str_vec_to_oss_str(std::vector<std::string> names) {
        int oss_name_len = max_strlen(names);
        std::stringstream oss_ss;
        for (const auto& name : names) {
            oss_ss << std::setw (oss_name_len) << name;
        }
        return oss_ss.str();
    }

    OssFixedInputs fixed_inputs;
    OssFixedOutputs fixed_outputs;
};

}



#endif
