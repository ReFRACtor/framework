#ifndef OSS_FORWARD_MODEL_H
#define OSS_FORWARD_MODEL_H
#include "forward_model.h"
#include "state_vector.h"
#include "rt_atmosphere.h"

#include <string>

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
        const std::string& Ch_sel_file, int Max_chans = 20000) : atmosphere(Atm),
          sel_file(Sel_file), sel_file_sz(Sel_file.length()),
          od_file(Od_file), od_file_sz(Od_file.length()),
          sol_file(Sol_file), sol_file_sz(Sol_file.length()),
          fix_file(Fix_file), fix_file_sz(Fix_file.length()),max_chans(Max_chans),
          ch_sel_file(Ch_sel_file), ch_sel_file_sz(Ch_sel_file.length()),
          max_chans(Max_chans) {} ;
    virtual ~OssForwardModel() {}
    virtual void setup_grid() {
      /*
      cppinitwrapper(nInMol, lenG, nameGas, nInJac, lenJ, nameJacob,
     sel_file, sel_file_sz,
     od_file, od_file_sz,
     sol_file, sol_file_sz,
     fix_file, fix_file_sz,
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
      char nameGas[110] = "H2O       CO2       O3        N2O       CO        CH4       O2        NH3       CCL4      F11       F12       ";
      int nInJac = 2;
      int lenJ = 2;
      char nameJacob[30] = "H2O       NH3       i";
      int nlevu = 65;
      int n_SfGrd = 501;
      float minExtCld = 999;
      int nchanOSS;
      float WvnOSS[max_chans];

      cppinitwrapper(nInMol, lenG, nameGas, nInJac, lenJ, nameJacob,
           sel_file.c_str(), sel_file_sz,
           od_file.c_str(), od_file_sz,
           sol_file.c_str(), sol_file_sz,
           fix_file.c_str(), fix_file_sz,
           ch_sel_file.c_str(), ch_sel_file_sz,
           nlevu, n_SfGrd, minExtCld,
           nchanOSS, WvnOSS, max_chans);

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

      /*
(lldb) p nlevu
(int) $0 = 65

(lldb) p ngas
(int) $1 = 11

tape5.nc/Pressure
(lldb) p Pin
(float *) $2 = 0x0000000100c530c0
(lldb) x/20f Pin
0x100c530c0: 1009.40991
0x100c530c4: 1000
0x100c530c8: 908.513977
0x100c530cc: 825.401978
0x100c530d0: 749.893005
0x100c530d4: 681.291015
0x100c530d8: 618.966003
0x100c530dc: 562.34198
0x100c530e0: 510.89801
0x100c530e4: 464.160004
0x100c530e8: 421.697998
0x100c530ec: 383.117004
0x100c530f0: 348.069
0x100c530f4: 316.22699
0x100c530f8: 287.298004
0x100c530fc: 261.015991
0x100c53100: 237.136993
0x100c53104: 215.444
0x100c53108: 195.735001
0x100c5310c: 177.828995

tape5.nc/Temperature
(lldb) p Temp
(float *) $3 = 0x0000000100c53200

tape5.nc/SkinTemperature
(lldb) p Tskin
(float) $4 = 310

tape5.nc/vmrGas
(lldb) p vmrGas
(float *) $5 = 0x000000010104f600

(lldb) p n_SfGrd
(int) $6 = 501

tape5.nc/Emissivity
(lldb) p emis
(float *) $7 = 0x000000010100d000

tape5.nc/Reflectivity
(lldb) p refl
(float *) $8 = 0x0000000101050200

(lldb) p scaleCld
(float) $9 = 0

(lldb) p presCld
(float) $10 = 0

(lldb) p nCld
(int) $11 = 2

(lldb) p extCld
(float *) $12 = 0x0000000100c531d0

tape5.nc/SurfaceGrid
(lldb) p SfGrd
(float *) $13 = 0x0000000101001800

(lldb) p gridCld
(float *) $14 = 0x0000000100c53080

(lldb) p ang
(float) $15 = 1.45646667

(lldb) p sunang
(float) $16 = 90

tape5.nc/Latitude
(lldb) p lat
(float) $17 = 45

(lldb) p altSfc
(float) $18 = 0.0000639999998

(lldb) p lambertian
(int) $19 = 1

(lldb) p nInJac
(int) $20 = 2

(lldb) p nchanOSS
(int) $21 = 3951

       */
      return Spectrum();
    }
    virtual void print(std::ostream& Os) const { Os << "OssForwardModel"; }
private:
    boost::shared_ptr<RtAtmosphere> atmosphere;
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

    const int max_chans;
};
}
#endif
