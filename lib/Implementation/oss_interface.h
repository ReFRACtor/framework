#ifndef OSS_INTERFACE_H
#define OSS_INTERFACE_H

#include <string>
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
    OssModifiedOutputs(blitz::Array<float, 1>& Y, blitz::Array<float, 1>& XkTemp, blitz::Array<float, 1>& XkTskin,
            blitz::Array<float, 1>& XkOutGas, blitz::Array<float, 1>& XkEm, blitz::Array<float, 1>& XkRf,
            blitz::Array<float, 1>& XkCldlnPres, blitz::Array<float, 1>& XkCldlnExt) :
            y(Y), xk_temp(XkTemp), xk_tskin(XkTskin), xk_out_gas(XkOutGas), xk_em(XkEm),
            xk_rf(XkRf), xk_cldln_pres(XkCldlnPres), xk_cldln_ext(XkCldlnExt){
    }

    const blitz::Array<float, 1>& getXkCldlnExt() const {
        return xk_cldln_ext;
    }

    const blitz::Array<float, 1>& getXkCldlnPres() const {
        return xk_cldln_pres;
    }

    const blitz::Array<float, 1>& getXkEm() const {
        return xk_em;
    }

    const blitz::Array<float, 1>& getXkOutGas() const {
        return xk_out_gas;
    }

    const blitz::Array<float, 1>& getXkRf() const {
        return xk_rf;
    }

    const blitz::Array<float, 1>& getXkTemp() const {
        return xk_temp;
    }

    const blitz::Array<float, 1>& getXkTskin() const {
        return xk_tskin;
    }

    const blitz::Array<float, 1>& getY() const {
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
    OssFixedOutputs(int NumChan, ArrayWithUnit<double, 1>& CenterWavenumber) :
            num_chan(NumChan), center_wavenumber(CenterWavenumber) {
    }

    const ArrayWithUnit<double, 1>& getCenterWavenumber() const {
        return center_wavenumber;
    }

    int getNumChan() const {
        return num_chan;
    }

private:
    int num_chan; ///< Number of channels available in OSS RTM
    ArrayWithUnit<double, 1> center_wavenumber; //< Center wavenumbers of channels

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

    const std::string& getChSelFile() const {
        return ch_sel_file;
    }

    const std::string& getFixFile() const {
        return fix_file;
    }

    const std::vector<std::string>& getGasJacobianNames() const {
        /* TODO: Translate to 1d C char array */
        return gas_jacobian_names;
    }

    const std::vector<std::string>& getGasNames() const {
        /* TODO: Translate to 1d char array */
        return gas_names;
    }

    int getMaxChans() const {
        return max_chans;
    }

    float getMinExtinctCld() const {
        return min_extinct_cld;
    }

    int getNumMolecules() const {
        return num_molecules;
    }

    int getNumSurfPoints() const {
        return num_surf_points;
    }

    int getNumVertLev() const {
        return num_vert_lev;
    }

    const std::string& getOdFile() const {
        return od_file;
    }

    const std::string& getSelFile() const {
        return sel_file;
    }

    const std::string& getSolFile() const {
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

    const blitz::Array<float, 1>& getCldGrid() const {
        return cld_grid;
    }

    void setCldGrid(const blitz::Array<float, 1>& cldGrid) {
        cld_grid = cldGrid;
    }

    const blitz::Array<float, 1>& getEmis() const {
        return emis;
    }

    void setEmis(const blitz::Array<float, 1>& emis) {
        this->emis = emis;
    }

    const blitz::Array<float, 1>& getExtCld() const {
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

    float getLat() const {
        return lat;
    }

    void setLat(float lat) {
        this->lat = lat;
    }

    float getObsZenAng() const {
        return obs_zen_ang;
    }

    void setObsZenAng(float obsZenAng) {
        obs_zen_ang = obsZenAng;
    }

    const blitz::Array<float, 1>& getPressure() const {
        return pressure;
    }

    void setPressure(const blitz::Array<float, 1>& pressure) {
        this->pressure = pressure;
    }

    float getPressureCld() const {
        return pressure_cld;
    }

    void setPressureCld(float pressureCld) {
        pressure_cld = pressureCld;
    }

    const blitz::Array<float, 1>& getRefl() const {
        return refl;
    }

    void setRefl(const blitz::Array<float, 1>& refl) {
        this->refl = refl;
    }

    float getScaleCld() const {
        return scale_cld;
    }

    void setScaleCld(float scaleCld) {
        scale_cld = scaleCld;
    }

    float getSkinTemp() const {
        return skin_temp;
    }

    void setSkinTemp(float skinTemp) {
        skin_temp = skinTemp;
    }

    float getSolZenAng() const {
        return sol_zen_ang;
    }

    void setSolZenAng(float solZenAng) {
        sol_zen_ang = solZenAng;
    }

    float getSurfAlt() const {
        return surf_alt;
    }

    void setSurfAlt(float surfAlt) {
        surf_alt = surfAlt;
    }

    const blitz::Array<float, 1>& getSurfGrid() const {
        return surf_grid;
    }

    void setSurfGrid(const blitz::Array<float, 1>& surfGrid) {
        surf_grid = surfGrid;
    }

    const blitz::Array<float, 1>& getTemp() const {
        return temp;
    }

    void setTemp(const blitz::Array<float, 1>& temp) {
        this->temp = temp;
    }

    const blitz::Array<float, 2>& getVmrGas() const {
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
        /* Allocate max_chans size arrays for OssFixedOutputs */
        /* Call cppinitwrapper */
        /* Resize arrays to returned number of channels */
        /* Set internal fixed_outputs */
    }

    OssModifiedOutputs run_fwd_model(OssModifiedInputs& ModifiedInputs) {
        /* Allocate memory for OssModifiedOutputs members */
        /* Call cppfwdwrapper */
        /* Cleanup / Conversion */
        return OssModifiedOutputs();
    }

private:
    OssFixedInputs fixed_inputs;
    OssFixedOutputs fixed_outputs;
};

}



#endif
