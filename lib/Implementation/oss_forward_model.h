#ifndef OSS_FORWARD_MODEL_H
#define OSS_FORWARD_MODEL_H
#include "forward_model.h"
#include "state_vector.h"
#include "rt_atmosphere.h"

namespace FullPhysics {
/****************************************************************//**
  This a forward model class that wraps the AER OSS Forward Model
*******************************************************************/
extern "C" {
  void cppinitwrapper(int &, int &, char *, int &, int &, char *,
                      const char *, int &,
                      const char *, int &,
                      const char *, int &,
                      const char *, int &,
                      const char *, int &,
                      int &, int &, float &,
                      int &, float *, const int &);
}

extern "C" {
  void cppfwdwrapper(int &, int &, float *, float *, float &, float *,
                     int &, float *, float *,
                     float &, float &, int &, float *,
                     float *, float *,
                     float &, float &, float &,
                     float &, int &, int &, int &,
                     float *, float *,
                     float *, float *, float *,
                     float *, float *, float *);
}

class OssForwardModel : public ForwardModel {
public:
    OssForwardModel(const boost::shared_ptr<RtAtmosphere>& Atm, const std::string& Sel_file,
        const std::string& Od_file, const std::string& Sol_file, const std::string& Fix_file,
        const std::string& Fix_file2, const std::string& Ch_sel_file) :
          atmosphere(Atm), sel_file(Sel_file), od_file(Od_file), sol_file(Sol_file),
          fix_file(Fix_file), fix_file2(Fix_file2), ch_sel_file(Ch_sel_file){} ;
    virtual ~OssForwardModel() {}
    virtual void setup_grid() {
      /*
      cppinitwrapper(nInMol, lenG, nameGas, nInJac, lenJ, nameJacob,
     sel_file, sel_file_sz,
     od_file, od_file_sz,
     sol_file, sol_file_sz,
     fix_file, fix_file_sz,
     fix_file2, fix_file2_sz,
     chSel_file, chSel_file_sz,
     nlevu, n_SfGrd, minExtCld,
     nchanOSS, WvnOSS, mxchan);
      where the arguments are defined as follows.
      . Inputs
      - nInMol = number of molecules in input profiles
      - nameGas = molecular gas names
      - nameJacob = molecular gas names for Jacobians
      - sel_file = file name of OSS nodes and weights
      - od_file = file name of optical property lookup table
      - sol_file = file name of solar radiance*
      - fix_file = file name of default profiles for variable gases
      - fix_file2 = file name of OSS data related to od_file (temporary)
      - chSelFile = file name of list(s) of channel subsets†
      - nlevu = number of vertical levels of state vector
      - nsf = number of surface grid points
      - minExtCld = threshold of extinction for including cloud in RT (km−1)

      . Outputs
      - nchanOSS = number of channels available in OSS RTM
      - cWvnOSS = center wavenumbers of channels
      */
      int nInMol = 11;
      int lenG = 10;
      char* nameGas[] = "H2O       CO2       O3        N2O       CO        CH4       O2        NH3       CCL4      F11       F12       ";
      int nInJac = 2;
      int lenJ = 2;
      char* nameJacob[] = "H2O       NH3       i";
      int nlevu = 65;
      int n_SfGrd = 501;
      float minExtCld = 999;
      int nchanOSS;
      float WvnOSS;
      const int mxchan;

    }
    virtual int num_channels() const { return 1; }
    virtual SpectralDomain spectral_domain(int Spec_index) const {
      return SpectralDomain();
    }
    virtual SpectralDomain::TypePreference spectral_domain_type_preference() const {
      return SpectralDomain::PREFER_WAVENUMBER;
    }
    virtual boost::shared_ptr<StateVector> state_vector() const 
    {
      return boost::shared_ptr<StateVector>();
    }

    virtual Spectrum radiance(int channel_index, bool skip_jacobian = false) const {
      /*
      cppfwdwrapper(nlevu, ngas, Pin, Temp, Tskin, vmrGas,
      n_SfGrd, emis, refl,
      scaleCld, presCld, nCld, extCld,
      SfGrd, gridCld,
      ang, sunang, lat,
      altSfc, lambertian, nInJac, nchan,
      y, xkTemp,
      xkTskin, xkOutGas, xkEm,
      xkRf, xkCldlnPres, xkCldlnExt);

      The I/Os are defined as follows.
      . Inputs
      - pin = atmospheric pressure profile (ordered from high to low) (mb)
      - temp = temperature profile (K)
      - tSkin = skin temperature (K)
      - vmrGas = volume mixing ratio of gases (unitless)
      - emis = surface emissivity (unitless)
      - refl = surface reflectivity (unitless)
      - scaleCld = scale of cloud log thickness
      - presCld = pressure at cloud center (peak of extinction profile) (mb)
      - extCld = cloud peak extinction (km−1)
      - SfcGrd = vector of surface grid point wavenumbers (cm−1)
      - gridCld = cloud grid point wavenumbers (cm−1)
      - ang = observation zenith angle (degrees)
      - sunang = solar zenith angle (degrees)
      - lat = latitude (degrees)
      - altSfc = surface altitude (m)
      - lambertian = flag for Lambertian (1) or (default) specular (0) surface

      . Outputs
      - y = radiance (W m−2 str−1 cm−1)
      - xkTemp = temperature profile Jacobians (W m−2 str−1 cm−1 K−1)
      - xkTskin = skin temperature Jacobians (W m−2 str−1 cm−1 K−1)
      - xkOutGas = gas log concentration Jacobians (W m−2 str−1 cm−1)
      - xkEm = emissivity Jacobians (W m−2 str−1 cm−1)
      - xkRF = reflectivity Jacobians (W m−2 str−1 cm−1)
      - xkCldlnPres = cloud center log pressure Jacobians  (W m−2 str−1 cm−1)
      - xkCldlnExt = cloud peak log extinction Jacobians  (W m−2 str−1 cm−1)
      */
      return Spectrum();
    }
    virtual void print(std::ostream& Os) const { Os << "OssForwardModel"; }
private:
    boost::shared_ptr<RtAtmosphere> atmosphere;
    std::string sel_file;
    std::string od_file;
    std::string sol_file;
    std::string fix_file;
    std::string fix_file2;
    std::string ch_sel_file;
};
}
#endif
