// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "fp_common.i"
%{
#include "ils_table.h"
%}
%base_import(ils_function)
%import "hdf_file.i"
%import "array_ad.i"
%import "auto_derivative.i"
%fp_shared_ptr(FullPhysics::IlsTableLinear);
%fp_shared_ptr(FullPhysics::IlsTableLog);

namespace FullPhysics {
class IlsTableLinear : public IlsFunction {
public:
  IlsTableLinear(const blitz::Array<double, 1>& Wavenumber, 
                 const blitz::Array<double, 2>& Delta_lambda, 
                 const blitz::Array<double, 2>& Response,
                 const std::string& Band_name, const std::string& Hdf_band_name,
                 bool Interpolate_wavenumber = false);
  IlsTableLinear(const HdfFile& Hdf_static_input, int Spec_index,
                 const std::string& Band_name, const std::string& Hdf_band_name,
                 const std::string& Hdf_group = "Instrument/ILS");
  virtual void ils
  (const AutoDerivative<double>& wn_center,
   const blitz::Array<double, 1>& wn, ArrayAd<double, 1>& OUTPUT) const;

  %python_attribute(band_name, std::string)
  %python_attribute(hdf_band_name, std::string)

  %python_attribute(wavenumber, blitz::Array<double, 1>)
  %python_attribute(delta_lambda, blitz::Array<double, 2>)
  %python_attribute(response, blitz::Array<double, 2>)

  void create_delta_lambda_to_response(const blitz::Array<double, 1>& Wavenumber, 
                                       const blitz::Array<double, 2>& Delta_lambda, 
                                       const blitz::Array<double, 2>& Response);
  %pickle_serialization();
};

class IlsTableLog : public IlsFunction {
public:
  IlsTableLog(const blitz::Array<double, 1>& Wavenumber, 
              const blitz::Array<double, 2>& Delta_lambda, 
              const blitz::Array<double, 2>& Response,
              const std::string& Band_name, const std::string& Hdf_band_name,
              bool Interpolate_wavenumber = false);
  IlsTableLog(const HdfFile& Hdf_static_input, int Spec_index,
              const std::string& Band_name, const std::string& Hdf_band_name,
              const std::string& Hdf_group = "Instrument/ILS");
  virtual void ils
  (const AutoDerivative<double>& wn_center,
   const blitz::Array<double, 1>& wn, ArrayAd<double, 1>& OUTPUT) const;

  %python_attribute(band_name, std::string)
  %python_attribute(hdf_band_name, std::string)

  %python_attribute(wavenumber, blitz::Array<double, 1>)
  %python_attribute(delta_lambda, blitz::Array<double, 2>)
  %python_attribute(response, blitz::Array<double, 2>)

  void create_delta_lambda_to_response(const blitz::Array<double, 1>& Wavenumber, 
                                       const blitz::Array<double, 2>& Delta_lambda, 
                                       const blitz::Array<double, 2>& Response);
  %pickle_serialization();
};

}
