#ifndef ILS_TABLE_H
#define ILS_TABLE_H
#include "ils_function.h"
#include "hdf_file.h"
#include "linear_interpolate.h"
#include "log_interpolate.h"

namespace FullPhysics {
/****************************************************************//**
  This class in a IlsFunction where we get the ILS response by using a
  table of measured values.

  The data is given to us in a fixed table of response functions,
  indexed by delta wavenumber (called delta_lambda in the HDF file),
  and for a small number of wavenumbers (called wavenumber in the HDF
  file). We always interpolate in the delta_lambda direction, but we
  may either pick the row of the table with the largest wavenumber <
  wn_center or we may interpolate between the two closest
  wavenumbers.

  In the old code, this was called KEY_ILS_TABLE and
  KEY_ILS_INTERP. In the Hdf file, this is "function_type" and is
  either "table" or "interpol".

  Note that "Interp" is a somewhat misleading name. We always
  interpolate in the delta lambda direction, the question is just if
  we also interpolate in the wavenumber direction.

  This one class handles both of the methods of interpolation,
  depending on the arguments passed to the function.
*******************************************************************/

template <class Lint>
class IlsTable : virtual public IlsFunction {
public:
//-----------------------------------------------------------------------
/// Constructor where we just supply the wavenumber, delta_lambda and
/// response values.
///
/// The arrays describe the ILS function. The array Wavenumber gives the
/// wavenumber of each pixel, Delta_lambda gives the difference
/// between the high resolution point we are looking at and the center
/// wavenumber, and response gives the detector response. 
///
/// The row gives values for a particular instrument wavenumber, for
/// GOSAT we typically have three rows for three different
/// wavenumbers, one on each edge and one at the middle. The columns
/// give a range of delta_lambda values.
//-----------------------------------------------------------------------

  IlsTable(const blitz::Array<double, 1>& Wavenumber, 
           const blitz::Array<double, 2>& Delta_lambda, 
           const blitz::Array<double, 2>& Response,
           const std::string& Band_name, const std::string& Hdf_band_name,
           bool Interpolate_wavenumber = false)
    : band_name_(Band_name), 
      hdf_band_name_(Hdf_band_name), 
      interpolate_wavenumber(Interpolate_wavenumber),
      from_hdf_file(false)
  { 
    create_delta_lambda_to_response(Wavenumber, Delta_lambda, Response); 
  }
  IlsTable(const HdfFile& Hdf_static_input, int Spec_index, 
           const std::string& Band_name, const std::string& Hdf_band_name,
           const std::string& Hdf_group = "Instrument/ILS");

  virtual ~IlsTable() {}
  virtual void ils
  (const AutoDerivative<double>& wn_center,
   const blitz::Array<double, 1>& wn,
   ArrayAd<double, 1>& res) const;

  virtual void print(std::ostream& Os) const;
  virtual std::string band_name() const {return band_name_;}
  virtual std::string hdf_band_name() const { return hdf_band_name_;}
  
  /// Accessors for retrieving arrays used to create table
  virtual blitz::Array<double, 1> wavenumber() const { return wavenumber_; }
  virtual blitz::Array<double, 2> delta_lambda() const { return delta_lambda_; }
  virtual blitz::Array<double, 2> response() const { return response_; }
  
  /// Creates the datastructures needed for setting up the ILS table,
  /// only necessary to rerun if modifying ILS table after object creation
  void create_delta_lambda_to_response(const blitz::Array<double, 1>& Wavenumber, 
                                       const blitz::Array<double, 2>& Delta_lambda, 
                                       const blitz::Array<double, 2>& Response);
private:
  void init_from_file(const HdfFile& Hdf_static_input);
  // Save some typing by giving these typedefs
  typedef ArrayAd<double, 1> arrad;
  typedef AutoDerivative<double> ad;
  typedef typename std::map<double, Lint>::const_iterator it;

  /// This is indexed by wavenumber, and returns a interpolator that
  /// goes from delta_lambda to the response.
  /// response arrays.
  std::map<double, Lint> delta_lambda_to_response;

  // Store the input arrays so we can re-build the ils map for testing
  // purposes, or for whatever reason you would need to do this mid program
  blitz::Array<double, 1> wavenumber_;
  blitz::Array<double, 2> delta_lambda_;
  blitz::Array<double, 2> response_;

  std::string band_name_, hdf_band_name_;
  bool interpolate_wavenumber;
  
  // Stash some information for print out with print
  bool from_hdf_file;
  std::string hdf_file_name;
  std::string hdf_group;
  IlsTable() {}
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
  template<class Archive>
  void save(Archive & ar, const unsigned int version) const;
  template<class Archive>
  void load(Archive & ar, const unsigned int version);
};

typedef IlsTable<LinearInterpolate<AutoDerivative<double>, AutoDerivative<double> > > IlsTableLinear;
typedef IlsTable<LinearLogInterpolate<AutoDerivative<double>, AutoDerivative<double> > > IlsTableLog;

}

FP_EXPORT_KEY(IlsTableLinear);
FP_EXPORT_KEY(IlsTableLog);
#endif
