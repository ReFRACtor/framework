#ifndef OSS_FORWARD_MODEL_H
#define OSS_FORWARD_MODEL_H
#include "forward_model.h"
#include "state_vector.h"
#include "rt_atmosphere.h"
#include "hdf_file.h"

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
          fix_file(Fix_file), fix_file_sz(Fix_file.length()),
          ch_sel_file(Ch_sel_file), ch_sel_file_sz(Ch_sel_file.length()),
          max_chans(Max_chans), nchanOSS(-1){} ;
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
      char nameGas[150] = "H2O       CO2       O3        N2O       CO        CH4       O2        NH3       CCL4      F11       F12       ";
      int nInJac = 2;
      int lenJ = 10;
      char nameJacob[30] = "H2O       NH3       i";
      int nlevu = 65;
      int n_SfGrd = 501;
      float minExtCld = 999;
      /* TODO: Make private data member */
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
      This method is primarily a wrapper for the following OSS FM function.

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
      using namespace blitz;

      const std::string test_fname("/export/menja_scr/refractor/OSS/run/tape5_nc4.nc");
      boost::shared_ptr<HdfFile> test_input_file(new HdfFile(test_fname));

      int nlevu = 65;

      int ngas = 11;

      Array<float, 1> Pin = test_input_file->read_field<float, 1>("/Pressure")(Range::all());

      Array<float, 1> Temp = test_input_file->read_field<float, 1>("/Temperature")(Range::all());

      float Tskin = 310.0;

      Array<float, 2> vmrGas = test_input_file->read_field<float, 2>("/vmrGas")(Range::all());

      Array<float, 1> emis = test_input_file->read_field<float, 1>("/Emissivity")(Range::all());

      Array<float, 1> refl = test_input_file->read_field<float, 1>("/Reflectivity")(Range::all());

      int n_SfGrd = 501;

      float scaleCld = 0.0;

      float presCld = 0.0;

      int nCld = 2;

      Array<float, 1> extCld(nCld);
      extCld = 0;

      Array<float, 1> SfGrd = test_input_file->read_field<float, 1>("/SurfaceGrid")(Range::all());

      Array<float, 1> gridCld(nCld);
      gridCld= 0;

      float ang = 1.45646667;

      float sunang = 90.0;

      float lat = 45.0;

      float altSfc = 0.0000639999998;

      int lambertian = 1;

      int nInJac = 2;

      /* Needed since nchanOSS argument to cppfwdwrapper isn't const reference though not modified */
      int nchanOSS_temp = nchanOSS;


      /* outputs */
      Array<float, 1> y(nchanOSS);
      y = 0;
      Array<float, 1> xkTemp(nchanOSS*nlevu);
      xkTemp = 0;
      Array<float, 1> xkTskin(nchanOSS);
      xkTskin = 0;
      Array<float, 1> xkOutGas(nInJac*nchanOSS*nlevu);
      xkOutGas = 0;
      Array<float, 1> xkEm(nchanOSS*n_SfGrd);
      xkEm = 0;
      Array<float, 1> xkRf(nchanOSS*n_SfGrd);
      xkRf = 0;
      Array<float, 1> xkCldlnPres(nchanOSS);
      xkCldlnPres = 0;
      Array<float, 1> xkCldlnExt(nchanOSS*nCld);
      xkCldlnExt = 0;

      cppfwdwrapper(nlevu, ngas, Pin.data(), Temp.data(), Tskin, vmrGas.data(),
      n_SfGrd, emis.data(), refl.data(),
      scaleCld, presCld, nCld, extCld.data(),
      SfGrd.data(), gridCld.data(),
      ang, sunang, lat,
      altSfc, lambertian, nInJac, nchanOSS_temp,
      y.data(), xkTemp.data(),
      xkTskin.data(), xkOutGas.data(), xkEm.data(),
      xkRf.data(), xkCldlnPres.data(), xkCldlnExt.data());

      SpectralDomain spec_domain;
      SpectralRange spec_range;
      return Spectrum(spec_domain, spec_range);
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

    int nchanOSS;
    int max_chans;
};
}
#endif
